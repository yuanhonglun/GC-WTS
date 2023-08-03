import pandas as pd
import numpy as np
import re
# import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
# import matplotlib.cm as cm
import math
from scipy.signal import find_peaks, peak_widths
from statistics import median
import pickle
import os
import pymzml
import gc
import operator


class DataAnalysis():

    ##################

    ##################
    # 本段将目标物RT±0.3 min中目标物的离子提取出来
    def find_nearest(self, array, value):
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or abs(value - array[idx - 1]) < abs(value - array[idx])):
            return idx - 1
        else:
            return idx

    ########
    # backup0215中第3部分原位置
    ########

    ##################

    ##################
    def calulate_lwma(self, x, weights, no_smooth_number):
        '''
            该函数的作用就是为了解决处理首位点的list_sum问题。由于local_df_padding是在local_df首尾拼接了两个为0的Dataframe，且该Dataframe的index是从0
            开始。思路就是取每次rolling得到的长度为n的series的index。如n=5，第一个rolling的index即为[0, 1, 899.554, 870.091, 870.628],切片拿到index[0]后
            可以推算出list_sum的长度
        '''

        if isinstance(x.index[0], str):
            list_sum = sum(weights[no_smooth_number - int(x.index[0]):])
            return np.dot(x, weights) / list_sum
        elif isinstance(x.index[-1], str):
            list_sum = sum(weights[int(x.index[-1]) + 1:])
            return np.dot(x, weights) / list_sum
        else:
            list_sum = sum(weights)
            return np.dot(x, weights) / list_sum

    ##################

    def validateTitle(self, title):
        rstr = r"[\/\\\:\*\?\"\<\>\|]"  # '/ \ : * ? " < > |'
        new_title = re.sub(rstr, "_", title)  # 替换为下划线
        return new_title

    ##################
    # 本段使用加权平均数平滑质谱图，使用了rolling、dot等快速方法可参考
    # 已添加振寰编写的边界值平滑方法，但未验证是否准确
    def lwma(self, local_df, n):
        n = int(n)
        smooth_df = pd.DataFrame(columns=local_df.columns.values, index=local_df.index.values)
        list1 = [i for i in range(1, int((n + 1) / 2) + 1)]
        weights = list1 + list(reversed(list1))[1:]
        no_smooth_number = int((n - 1) / 2)
        tmp_df = pd.DataFrame(np.zeros((no_smooth_number, local_df.shape[1])), columns=list(local_df.columns))
        tmp_df = tmp_df.set_axis(list(map(str, list(tmp_df.index))))
        # 将0矩阵axis变成字符串，判断是否为字符串
        local_df_padding = pd.concat([tmp_df, local_df, tmp_df])
        for ion in local_df.columns.values:
            smooth = local_df_padding[ion].rolling(window=n, min_periods=n, center=True).apply(self.calulate_lwma,
                                                                                               args=(
                                                                                                   weights,
                                                                                                   no_smooth_number)).to_list()
            # 此处window = 5, center = True时，即取某点及其前后各两点计算，min_periods = 1时，第一个点取其自己与后两点
            # (n = 5)的情况下
            smooth_df[ion] = smooth[no_smooth_number:-no_smooth_number]

        smooth_df = smooth_df.dropna(axis=0, how='all')

        return smooth_df

    ##################

    ##################
    # 本段计算出每个点的fd、sd和ad(amplitude difference)，部分点无法参与计算，暂时去掉
    # 计算出fd等数据后，顺带计算出了ff、sf、af等噪音阈值，结果经excel确认是准确的
    # np.median()取的中位数逻辑如下：如果序列是奇数，则取中位数，如果序列是偶数，则取中间两个数的平均值
    def derivative(self, smooth_df, m):
        noise_list = []
        noise_col = smooth_df.columns.values
        for ion in smooth_df.columns.values:
            fd = smooth_df[ion].rolling(window=5, min_periods=5, center=True).apply(
                lambda x: ((-2 * x[0]) - x[1] + x[3] + (2 * x[4]))
                          / 10, raw=True).to_list()
            abs_fd = [abs(i) for i in fd]
            sd = smooth_df[ion].rolling(window=5, min_periods=5, center=True).apply(
                lambda x: ((2 * x[0]) - x[1] - (2 * x[2]) - x[3] + (2 * x[4]))
                          / 7, raw=True).to_list()
            neg_sd = [i * -1 for i in sd]
            ad = smooth_df[ion].rolling(window=2, min_periods=2).apply(lambda x: x[1] - x[0], raw=True).to_list()
            abs_ad = smooth_df[ion].rolling(window=2, min_periods=2).apply(lambda x: abs(x[1] - x[0]),
                                                                           raw=True).to_list()

            abs_fd_1 = [x for x in abs_fd if np.isnan(x) == False]
            neg_sd_1 = [x for x in neg_sd if np.isnan(x) == False]
            abs_ad_1 = [x for x in abs_ad if np.isnan(x) == False]

            ##################
            # 下面方法为取序列中小于最大值5%的数，取其中位数作为噪音阈值，经验证下列代码才是MSDIAL思路
            ff = np.median(np.asarray([i for i in abs_fd_1 if 0 < i < max(abs_fd_1) * 0.1]))
            sf = np.median(np.asarray([i for i in neg_sd_1 if i < 0 < abs(i) and abs(i) < max(neg_sd_1) * 0.1]))
            af = np.median(np.asarray([i for i in abs_ad_1 if 0 < i < max(abs_ad_1) * 0.1]))

            ff = max(ff, m)
            if af < m:
                af = m
            # 若ff和af接近0，则取m，参数可改

            # noise_df.loc["ff", ion] = ff
            # noise_df.loc["sf", ion] = sf
            # noise_df.loc["af", ion] = af
            noise_tmp = [ff, sf, af]
            noise_list.append(noise_tmp)
            smooth_df = pd.concat([smooth_df, pd.DataFrame({str(ion) + "_fd": fd}, index=smooth_df.index)], axis=1)
            smooth_df = pd.concat([smooth_df, pd.DataFrame({str(ion) + "_abs_fd": abs_fd}, index=smooth_df.index)],
                                  axis=1)
            smooth_df = pd.concat([smooth_df, pd.DataFrame({str(ion) + "_sd": sd}, index=smooth_df.index)], axis=1)
            smooth_df = pd.concat([smooth_df, pd.DataFrame({str(ion) + "_neg_sd": neg_sd}, index=smooth_df.index)],
                                  axis=1)
            smooth_df = pd.concat([smooth_df, pd.DataFrame({str(ion) + "_ad": ad}, index=smooth_df.index)], axis=1)
            smooth_df = pd.concat([smooth_df, pd.DataFrame({str(ion) + "_abs_ad": abs_ad}, index=smooth_df.index)],
                                  axis=1)

        smooth_df = smooth_df.dropna(axis=0, how='any')
        noise_df = pd.DataFrame(noise_list, columns=["ff", "sf", "af"], index=noise_col).T

        return smooth_df, noise_df

    ##################

    def find_peak(self, ion_list, smooth_df, noise_df, filter_factor):
        peak_dic = {}
        con_peak_dic = {}
        NF = self.pre_de_redundancy_peak(smooth_df, ion_list)

        # 本段采用了老写法，最好用rolling
        for ion in ion_list:
            # print("ion = ", ion)
            peak_df = pd.DataFrame(columns=["left", "apex", "right", "SV"])
            ff = noise_df.loc["ff", ion]
            sf = noise_df.loc["sf", ion]
            # af = noise_df.loc["af", ion]
            f_list = smooth_df[str(ion) + "_fd"].values.tolist()
            # abs_f_list = smooth_df[str(ion) + "_abs_fd"].values.tolist()
            s_list = smooth_df[str(ion) + "_sd"].values.tolist()
            a_list = smooth_df[str(ion) + "_ad"].values.tolist()
            # abs_a_list = smooth_df[str(ion) + "_abs_ad"].values.tolist()
            rt_list = smooth_df.index.values  # 改
            i_list = smooth_df[ion].values.tolist()

            ###########################

            m = 0
            total_list = []
            index_list = []
            aa = 1
            while m <= len(f_list) - 3:
                peak_list = []
                # len比最大index大1，所以减3
                if aa != m:
                    aa = m
                    if f_list[m] > ff * filter_factor and f_list[m + 1] > ff * filter_factor and i_list[m] > 0 and i_list[
                        m + 1] > 0:
                        # print("m = ", m)
                        w = m - 2
                        if w < 0:
                            w = 0
                        tmp_list = i_list[w: m + 3]
                        # print("tmp_list = ", tmp_list)
                        u = tmp_list.index(i_list[m])
                        # print("u = ", u)
                        q = tmp_list.index(min(tmp_list))
                        # print("q = ", q)
                        p = m - (u - q)
                        # print("p = ", p)
                        if i_list[p] == 0:
                            p = m
                            left = rt_list[p]
                        else:
                            left = rt_list[p]
                            # 加上这个条件，如果5点法找最小点找到了0值，那么就用原左点m
                        peak_list.append(left)
                        # peak_list.append(rt_list[m])
                        # print("---------左", left)

                        for x in range(p + 1, len(f_list) - 3):
                            # print("p+1 = ", p + 1)
                            # print("len(f_list) - 3 = ", len(f_list) - 3)

                            # 冲顶峰峰顶的判断条件
                            if f_list[x] == f_list[x + 1] == f_list[x + 2] == 0 and f_list[x + 3] != 0:

                                peak_list.append(rt_list[x])
                                index_list.append("peak_" + str(rt_list[x + 3]))
                                x += 4
                                for y in range(x, len(f_list) - 3):
                                    if f_list[y] > -ff * 10 and f_list[y + 1] > -ff * 10:
                                        w_1 = y + 3
                                        if w_1 > len(rt_list) + 1:
                                            w_1 = len(rt_list) + 1
                                        tmp_list = i_list[y - 2: w_1]
                                        u_1 = tmp_list.index(i_list[y])
                                        q_1 = tmp_list.index(min(tmp_list))
                                        p_1 = y - (u_1 - q_1)

                                        if i_list[p_1] == 0:
                                            right = rt_list[y]
                                            m = y
                                        elif rt_list[p_1] <= rt_list[x]:
                                            right = rt_list[y]
                                            m = y
                                        else:
                                            right = rt_list[p_1]
                                            m = p_1
                                        peak_list.append(right)
                                        # print("---------右", right)

                                        break

                                    elif i_list[y] < i_list[x] * 0.05:
                                        if rt_list[y] <= rt_list[x]:
                                            peak_list.append(rt_list[y + 1])
                                            m = y + 1
                                            # print("---------右", rt_list[y+1])
                                        else:
                                            peak_list.append(rt_list[y])
                                            m = y
                                            # print("---------右", rt_list[y])
                                        break
                                break
                            elif x > len(f_list) - 6:
                                peak_list = []
                                m = x
                                break

                            else:
                                if f_list[x - 1] > 0 and s_list[x] < -1 * sf:
                                    if f_list[x] < 0 or f_list[x + 1] < 0:

                                        # 找到顶点后，找其左右2点及本身，取最大值
                                        start_index = max(0, x - 2)
                                        end_index = min(len(i_list), x + 3)
                                        five_elements = i_list[start_index:end_index]
                                        max_index = five_elements.index(max(five_elements))
                                        n = start_index + max_index
                                        # print(rt_list[n])
                                        peak_list.append(rt_list[n])
                                        index_list.append("peak_" + str(rt_list[n]))

                                        for y in range(n + 1, len(f_list) - 3):
                                            if f_list[y] > -ff * filter_factor and f_list[y + 1] > -ff * filter_factor:
                                                w_1 = y + 3
                                                if w_1 > len(rt_list) + 1:
                                                    w_1 = len(rt_list) + 1
                                                tmp_list = i_list[y - 2: w_1]
                                                # print("tmp_list = ", tmp_list)
                                                u_1 = tmp_list.index(i_list[y])
                                                # print("u_1 = ", u_1)
                                                q_1 = tmp_list.index(min(tmp_list))
                                                # print("q_1 = ", q_1)
                                                p_1 = y - (u_1 - q_1)
                                                # print("p_1 = ", p_1)
                                                if i_list[p_1] == 0:
                                                    right = rt_list[y]
                                                    m = y
                                                elif rt_list[p_1] <= rt_list[n]:
                                                    right = rt_list[y]
                                                    m = y
                                                else:
                                                    right = rt_list[p_1]
                                                    m = p_1
                                                peak_list.append(right)

                                                # print("---------右282", right)
                                                break

                                            elif i_list[y] < i_list[n] * 0.05:
                                                if rt_list[y] <= rt_list[n]:
                                                    peak_list.append(rt_list[y + 1])
                                                    m = y + 1
                                                    # print("---------右289", rt_list[y+1])
                                                else:
                                                    peak_list.append(rt_list[y])
                                                    m = y
                                                    # print("---------右293", rt_list[y])
                                                break
                                            elif y == len(f_list) - 4:
                                                peak_list.append(rt_list[y])
                                                m = y
                                                break

                                        break
                        if peak_list != []:
                            total_list.append(peak_list)
                    else:
                        m += 1
                else:
                    print('find_peak_error')
                    break

            print(total_list)

            # 补充一个index_list作为行索引
            # 之后dataframe要添加一个空的SV列[[]]
            peak_df = pd.DataFrame(total_list, columns=["left", "apex", "right"], index=index_list)
            # peak_df.to_excel('./results/peak/{}.xlsx'.format(str(ion) + "_raw_peak"))
            peak_df = self.trailing_peak_filtration(peak_df, smooth_df, ion)
            # peak_df.to_excel('./results/peak/{}.xlsx'.format(str(ion) + "_filtered_peak"))
            # ##################画图#########################
            # for peak in peak_df.index.values:
            #     # print("debug_plot_ion = ", ion)
            #     apex_rt = peak_df.loc[peak, "apex"]
            #     apex_i = smooth_df.loc[apex_rt, ion]
            #     left_rt = peak_df.loc[peak, "left"]
            #     left_i = smooth_df.loc[left_rt, ion]
            #     right_rt = peak_df.loc[peak, "right"]
            #     right_i = smooth_df.loc[right_rt, ion]
            #     # print("debug_axe", apex_rt, apex_i, left_rt, left_i, right_rt, right_i)
            #     x = np.asarray(rt_list)
            #     y = np.asarray(i_list)
            #     ax = plt.axes()
            #     ax.plot(x, y, 'o')
            #     ax.set(xlabel="RT", ylabel="Intensity", title=ion)
            #     plt.annotate("*", (left_rt, left_i), xytext=(left_rt, left_i))
            #     plt.annotate("*", (apex_rt, apex_i), xytext=(apex_rt, apex_i))
            #     plt.annotate("*", (right_rt, right_i), xytext=(right_rt, right_i))
            #     plt.ylim(min(left_i, right_i) * 0.8, apex_i * 1.5)
            #     plt.xlim(left_rt - 5, right_rt + 5)
            #     plt.savefig('./results/peak/{}.png'.format(str(ion) + "_" + str(apex_rt)))
            #     plt.show()
            # ##################画图#########################
            ###基线矫正
            fb = 'smooth_df_final' in dir()
            if fb == False:
                from copy import copy  # 改
                smooth_df_final = copy(smooth_df.applymap(lambda x: x if x >= 0 else 0))
            # smooth_df_final.to_excel('./results/smooth_df_final_before_ba.xlsx')
            smooth_df_final = self.baseline_correction(peak_df, smooth_df, ion, smooth_df_final)
            # smooth_df_final.to_excel('./results/smooth_df_final_after_ba.xlsx')
            peak_df = self.de_redundancy_peak(smooth_df_final, smooth_df, peak_df, NF, ion)

            #########################################################
            ###3点抛物线校正
            rt_list = rt_list.tolist()
            for peak in peak_df.index.values:
                apex_rt = float(peak_df.loc[peak, "apex"])
                apex_i = float(smooth_df_final.loc[apex_rt, ion])
                apex_rt_index = rt_list.index(apex_rt)

                apex_left_rt = float(rt_list[apex_rt_index - 1])
                apex_left_i = float(smooth_df_final.loc[apex_left_rt, ion])
                apex_right_rt = float(rt_list[apex_rt_index + 1])
                apex_right_i = float(smooth_df_final.loc[apex_right_rt, ion])

                # 顶点及左右两点，共三点拟合抛物线，抛物线一阶导数为0的点为新顶点，计算出新的三点
                new_apex_rt, new_apex_i = self.second_order_fit_X(
                    apex_left_rt, apex_left_i, apex_rt, apex_i, apex_right_rt, apex_right_i)

                # 拟合后现将原来点的响应替换为nan
                # smooth_df_final.loc[apex_left_rt, ion] = np.nan
                smooth_df_final.loc[apex_rt, ion] = np.nan
                # smooth_df_final.loc[apex_right_rt, ion] = np.nan

                # 再填入拟合后的点，会使RT轴多出新点，其他离子在这些点的响应均为nan，但后续计算dropna会去掉
                # smooth_df_final.loc[new_apex_left_rt, ion] = new_apex_left_i
                smooth_df_final.loc[new_apex_rt, ion] = new_apex_i
                # smooth_df_final.loc[new_apex_rt, ion] = apex_i
                # smooth_df_final.loc[new_apex_right_rt, ion] = new_apex_right_i

                smooth_df_final.sort_index(inplace=True)

                # 最后将peak_df里的apex替换掉
                peak_df.loc[peak, "apex"] = new_apex_rt

                left_rt = peak_df.loc[peak, "left"]
                left_i = smooth_df_final.loc[left_rt, ion]
                right_rt = peak_df.loc[peak, "right"]
                right_i = smooth_df_final.loc[right_rt, ion]
                x_df = smooth_df_final.loc[left_rt: right_rt, ion]
                x_df = x_df.dropna(axis=0, how="any")
                # print("x_df = ", x_df)
                y = x_df.values  # 改
                # ax = plt.axes()
                # plt.plot(x_df.index, y, 'o')
                # ax.set(xlabel="RT", ylabel="Intensity", title=ion)
                # plt.annotate("*", (left_rt, left_i), xytext=(left_rt, left_i))
                # plt.annotate("*", (new_apex_rt, new_apex_i), xytext=(new_apex_rt, new_apex_i))
                # plt.annotate("*", (right_rt, right_i), xytext=(right_rt, right_i))
                # plt.ylim(min(left_i, right_i) * 0.8, apex_i * 1.5)
                # plt.xlim(left_rt - 5, right_rt + 5)
                # #plt.savefig('./results/peak/{}.png'.format(str(ion) + "_" + str(apex_rt) + "_adjusted_apex_RT"))
                # plt.close()
                del x_df, apex_left_rt, apex_right_rt
                gc.collect()
            # smooth_df_final.to_excel('./results/smooth_df_final_after_3_point.xlsx')
            # # peak_df.to_excel('./results/peak/{}.xlsx'.format(str(ion) + "_baseline_adjusted_peak"))
            ##################

            ##################
            # 1221卷积峰判断及提取算法：由于基线校正后，avl, avr均很小，因此舍弃该算法
            # 同时fd是粗放计算方法，导致肉眼看到fd为负时，计算出的fd为正，因此暂认为ideal_slope判断不出需要解卷积的峰, 及其local_apex
            # 现采用ad判断某峰某一边是否需要解卷积，左边:连续3点ad为正负负，且其ad均>af*5,那么第1点为左边界的loacl_apex
            # 右边：连续3点ad为正正负，且其ad均>af*5,那么第2点为右边界的loacl_apex
            # 连续几点判断，以及af乘以的系数，均可根据实际结果返回调整

            # 把该离子当前响应为NA的时间点先去掉
            temp_df = smooth_df_final[[ion]]  # 改
            # print("temp_df = ", temp_df)
            temp_df = temp_df.dropna(axis=0, how="any")
            i_list_1 = temp_df[ion].values.tolist()
            rt_list_1 = temp_df.index.values.tolist()

            con_peak_df = pd.DataFrame(columns=["intensity", "rt", "width", "raw_left", "raw_right"])
            for peak in peak_df.index.values:

                left = peak_df.loc[peak, "left"]
                apex = peak_df.loc[peak, "apex"]
                right = peak_df.loc[peak, "right"]
                left_index = rt_list_1.index(left)
                apex_index = rt_list_1.index(apex)
                right_index = rt_list_1.index(right)

                n = left_index + 1
                while left_index + 1 <= n <= apex_index - 3:
                    # todo: 确认一下找卷积峰条件
                    if a_list[n] > 0 and a_list[n + 1] > 0 and a_list[n + 2] < 0 and s_list[n + 1] < sf and i_list_1[
                        n + 1] >= i_list_1[apex_index] * 0.1:
                        if rt_list[n + 1] not in con_peak_df["rt"].values.tolist():
                            con_peak_df.loc["peak_" + str(rt_list[n + 1]), "intensity"] = i_list_1[n + 1]
                            con_peak_df.loc["peak_" + str(rt_list[n + 1]), "rt"] = rt_list[n + 1]
                            con_peak_df.loc["peak_" + str(rt_list[n + 1]), "width"] = (
                                                                                              right - left) / 2
                            con_peak_df.loc["peak_" + str(rt_list[n + 1]), "raw_left"] = left
                            con_peak_df.loc["peak_" + str(rt_list[n + 1]), "raw_right"] = right
                        if apex not in con_peak_df["rt"].values.tolist():
                            con_peak_df.loc[peak, "intensity"] = i_list_1[apex_index]
                            con_peak_df.loc[peak, "rt"] = apex
                            con_peak_df.loc[peak, "width"] = (right - left) / 2
                            con_peak_df.loc[peak, "raw_left"] = left
                            con_peak_df.loc[peak, "raw_right"] = right
                            if apex in peak_df["apex"].values.tolist():
                                peak_df = peak_df.drop(peak, axis=0)
                        n = n + 1
                    else:
                        n = n + 1

                m = apex_index + 1
                while apex_index + 1 <= m <= right_index - 3:
                    if a_list[m] < 0 and a_list[m + 1] < 0 and a_list[m + 2] > 0 and s_list[m + 2] < sf and i_list_1[
                        m + 2] >= i_list_1[apex_index] * 0.1:
                        if rt_list[m + 2] not in con_peak_df["rt"].values.tolist():
                            con_peak_df.loc["peak_" + str(rt_list[m + 2]), "intensity"] = i_list_1[m + 2]
                            con_peak_df.loc["peak_" + str(rt_list[m + 2]), "rt"] = rt_list[m + 2]
                            con_peak_df.loc["peak_" + str(rt_list[m + 2]), "width"] = (
                                                                                              right - left) / 2
                            con_peak_df.loc["peak_" + str(rt_list[m + 2]), "raw_left"] = left
                            con_peak_df.loc["peak_" + str(rt_list[m + 2]), "raw_right"] = right
                        if apex not in con_peak_df["rt"].values.tolist():
                            con_peak_df.loc[peak, "intensity"] = i_list_1[apex_index]
                            con_peak_df.loc[peak, "rt"] = apex
                            con_peak_df.loc[peak, "width"] = (right - left) / 2
                            con_peak_df.loc[peak, "raw_left"] = left
                            con_peak_df.loc[peak, "raw_right"] = right
                            if apex in peak_df["apex"].values.tolist():
                                peak_df = peak_df.drop(peak, axis=0)
                        m = m + 1
                    else:
                        m = m + 1

            con_peak_df.sort_values(by="rt", inplace=True, ascending=True)
            if not con_peak_df.empty:
                drop_con_peak_df_list = []
                for row in con_peak_df.itertuples():
                    if getattr(row, 'raw_left') > getattr(row, 'rt') or getattr(row, 'rt') > getattr(row,
                                                                                                     'raw_right') or getattr(
                        row, 'intensity') < 0:
                        drop_con_peak_df_list.append(row[0])
                con_peak_df = con_peak_df.drop(drop_con_peak_df_list)

            # con_peak_df.to_excel('./results/decon/{}.xlsx'.format(str(ion) + "_con_peak"))

            con_peak_dic[ion] = con_peak_df

            peak_dic[ion] = peak_df

            ##################

        return peak_dic, con_peak_dic, smooth_df_final

    def trailing_peak_filtration(self, peak_df, smooth_df, ion):
        """
        删除拖尾峰和识别左右错误的峰
        """
        for peak in peak_df.index.values:
            apex_rt = peak_df.loc[peak, "apex"]
            apex_i = smooth_df.loc[apex_rt, ion]
            left_rt = peak_df.loc[peak, "left"]
            left_i = smooth_df.loc[left_rt, ion]
            right_rt = peak_df.loc[peak, "right"]
            right_i = smooth_df.loc[right_rt, ion]
            if left_rt < apex_rt < right_rt:
                if right_i > 0.85 * apex_i or left_i > 0.85 * apex_i:
                    peak_df = peak_df.drop(peak, axis=0)
                else:
                    pass
            else:
                peak_df = peak_df.drop(peak, axis=0)

        return peak_df

    ##################

    ##找精确apex RT并拟合顶点及两边点的RT和intensity
    def second_order_fit_X(self, x1, y1, x2, y2, x3, y3):
        if x1 < x2 < x3 and y2 > y1 and y2 > y3:

            X = np.array([x1, x2, x3])
            Y = np.array([y1, y2, y3])
            coef = np.polyfit(X, Y, 2)
            x2_new = (-coef[1]) / (2 * coef[0])

            if x1 < x2_new < x3:
                y2_new = coef[0] * (x2_new ** 2) + coef[1] * x2_new + coef[2]
                if y2_new > y1 and y2_new > y3 and y2_new <= 1.5 * y2:
                    # 部分峰型奇怪的峰，拟合后y2_new会过大，这种情况返回原值
                    # 判断系数是1.5，此处修改
                    y2_new = y2_new
                else:
                    y2_new = y2
            else:
                x2_new = x2
                y2_new = y2
            # delta_x = x2_new - x2
            # x1_new = x1 + delta_x
            # x3_new = x3 + delta_x
            # y1_new = coef[0] * (x1_new ** 2) + coef[1] * x1_new + coef[2]
            # y3_new = coef[0] * (x3_new ** 2) + coef[1] * x3_new + coef[2]

            # return x1_new, y1_new, x2_new, y2_new, x3_new, y3_new
            return x2_new, y2_new

        else:

            return x2, y2

    ##################

    ##################
    def Gaussian(self, x, *params):
        num_func = int(len(params) / 3)
        y_list = []
        for i in range(num_func):
            y = np.zeros_like(x)
            param_range = list(range(3 * i, 3 * (i + 1), 1))
            amp = abs(params[int(param_range[0])])
            ctr = abs(params[int(param_range[1])])
            wid = abs(params[int(param_range[2])])
            y = y + amp * np.exp(-((x - ctr) / wid) ** 2)
            y_list.append(y)
        y_sum = np.zeros_like(x)
        for i in y_list:
            y_sum = y_sum + i
        y_sum = y_sum + params[-1]

        return y_sum

    def fit_plot(self, x, *params):
        num_func = int(len(params) / 3)
        y_list = []
        for i in range(num_func):
            y = np.zeros_like(x)
            param_range = list(range(3 * i, 3 * (i + 1), 1))
            amp = abs(params[int(param_range[0])])
            ctr = abs(params[int(param_range[1])])
            wid = abs(params[int(param_range[2])])
            y = y + amp * np.exp(-((x - ctr) / wid) ** 2) + params[-1]
            y_list.append(y)

        return y_list

    def decon(self, ion_list, smooth_df_final, con_peak_dic):
        # 3点校正后，高斯拟合出来的峰，其apex_RT与原峰（已3点校正）基本一致，不需要再次3点校正

        decon_peak_dic = {}
        decon_data_df = pd.DataFrame(index=smooth_df_final.index)
        for ion in ion_list:
            if con_peak_dic[ion].shape[0] == 0:
                decon_peak_dic[ion] = pd.DataFrame()
            else:
                temp_df = smooth_df_final.loc[:, ion]  # 改
                temp_df = temp_df.dropna(axis=0, how="any")
                i_list = temp_df.values.tolist()
                rt_list = temp_df.index.values.tolist()
                raw_left_list = con_peak_dic[ion]["raw_left"].values.tolist()
                new_raw_left_list = list(set(raw_left_list))
                new_raw_left_list.sort(key=raw_left_list.index)
                decon_peak_dic[ion] = pd.DataFrame()

                for left_value in new_raw_left_list:
                    # print("debug_left = ", left_value)
                    para_list = []
                    # wid_g = []
                    bound_min = []
                    bound_max = []
                    for peak in con_peak_dic[ion].index.values:
                        if con_peak_dic[ion].loc[peak, "raw_left"] == left_value:
                            para_list.append(con_peak_dic[ion].loc[peak, "intensity"])
                            para_list.append(con_peak_dic[ion].loc[peak, "rt"])
                            para_list.append(con_peak_dic[ion].loc[peak, "width"])
                            left = con_peak_dic[ion].loc[peak, "raw_left"]
                            right = con_peak_dic[ion].loc[peak, "raw_right"]
                            left_index = rt_list.index(left)
                            right_index = rt_list.index(right)
                            x = rt_list[left_index: right_index + 1]
                            y = i_list[left_index: right_index + 1]
                            bound_min.append(con_peak_dic[ion].loc[peak, "intensity"] * 0.9)
                            bound_max.append(con_peak_dic[ion].loc[peak, "intensity"] * 1.2)
                            bound_min.append(con_peak_dic[ion].loc[peak, "rt"] - 1)
                            bound_max.append(con_peak_dic[ion].loc[peak, "rt"] + 1)
                            bound_min.append(0)
                            bound_max.append(con_peak_dic[ion].loc[peak, "width"])

                    y_0 = []
                    for i in y:
                        if i != 0:
                            y_0.append(i)
                    # y_0.append(0)
                    before_half_y = y_0[:(len(y_0) // 2)]
                    after_half_y = y_0[(len(y_0) // 2):]
                    background_min = max(min(y_0), 0)
                    background_max = max(min(before_half_y), min(after_half_y), 0.01)  # 因为上下限不能相同所以使用0.1避免重复。--jfyd
                    # print(background_min, '*********', background_max)
                    if background_min < background_max:
                        para_list.append(background_min)
                        bound_min.append(background_min)
                        bound_max.append(background_max)
                    else:
                        para_list.append(0)
                        bound_min.append(0)
                        bound_max.append(1000000000)

                    # print("debug_b_min = ", bound_min)
                    # print("debug_b_max = ", bound_max)
                    popt, pcov = curve_fit(self.Gaussian, x, y, bounds=(bound_min, bound_max), p0=para_list,
                                           maxfev=100000)
                    # 此处要添加，如果拟合次数超过10000，仍无法拟合，那么返回原始数据作为y_list
                    fit = self.Gaussian(x, *popt)
                    y_list = self.fit_plot(x, *popt)
                    # print("debug_y_list = ", y_list)
                    # ax = plt.axes()
                    # ax.set(xlabel="RT", ylabel="Intensity", title=ion)
                    # plt.plot(x, y, lw=2)
                    # plt.plot(x, fit, ls='-', c='black', lw=0)
                    baseline = np.zeros_like(x) + popt[-1]
                    # print('baseline', baseline)
                    # for n, i in enumerate(y_list):
                    #     plt.fill_between(x, i, baseline, facecolor=cm.rainbow(n / len(y_list)), alpha=0.6)
                    # #plt.savefig('./results/decon/{}.jpg'.format(str(ion) + "_" + str(left) + "_group_decon"))
                    # plt.close()

                    for peak in y_list:
                        df = pd.DataFrame()

                        peak_1 = peak.tolist()
                        if max(peak_1) - min(peak_1) > 10:
                            peaks, _ = find_peaks(peak_1, height=0)
                            if len(peaks) > 0:
                                wid = peak_widths(peak_1, peaks, rel_height=0.999)
                                new_apex = x[peaks[0]]
                                new_left = x[round(wid[2][0])]
                                new_right = x[round(wid[3][0])]
                                peak_name = str(ion) + "_peak_" + str(new_apex)
                                decon_peak_dic[ion].loc[peak_name, "left"] = new_left
                                decon_peak_dic[ion].loc[peak_name, "apex"] = new_apex
                                decon_peak_dic[ion].loc[peak_name, "right"] = new_right

                                # 使用pandas.DataFrame.update()函数将dataframe_2中a列的值更新到dataframe_1中对应位置的a列上
                                # dataframe_2的index包含dataframe_1
                                df = pd.DataFrame(peak, index=x, columns=[peak_name])
                                # decon_data_df.update(df)
                                decon_data_df = decon_data_df.combine_first(df)

                    # decon_peak_dic[ion].to_excel('./results/decon/{}.xlsx'.format(str(ion) + "_decon_peak"))

        return decon_peak_dic, decon_data_df

    ##################

    ##################
    def sharpness(self, peak_df, data_df):
        data_df["number_index"] = [x for x in range(len(data_df))]
        # print("peak_df = ", peak)
        SV_list = []
        for index1, row1 in peak_df.iterrows():
            apex = data_df.loc[row1["apex"]]
            right = data_df.loc[row1["right"]]
            left = data_df.loc[row1["left"]]
            ion = data_df.columns[0]
            left_list = []
            right_list = []
            for index2, row2 in data_df.iterrows():
                if left["number_index"] <= row2['number_index'] < apex["number_index"]:
                    left_list.append(
                        (apex[ion] - row2[ion]) / (apex['number_index'] - row2['number_index']) / apex[ion] ** 0.5)
                elif apex["number_index"] < row2['number_index'] <= right["number_index"]:
                    right_list.append(
                        (apex[ion] - row2[ion]) / (row2['number_index'] - apex['number_index']) / apex[ion] ** 0.5)

            if left_list and right_list:
                sharpness = (max(left_list) + max(right_list)) / 2
            elif left_list:
                # print(peak_df)
                # print(apex)
                # print("该点右list为空")
                sharpness = max(left_list)
            elif right_list:
                # print(peak_df)
                # print(apex)
                # print("该点左list为空")
                sharpness = max(right_list)
            else:
                # print(peak_df)
                # print(apex)
                # print("该点左右list为空")
                sharpness = 0

            SV_list.append(sharpness)

        peak_df["SV"] = SV_list

        return peak_df

    def calculate_sv(self, smooth_df_final, peak_dic, decon_data_df, decon_peak_dic):
        keys_to_delete = []

        for key in peak_dic:
            if peak_dic[key].empty:
                keys_to_delete.append(key)

        for key in keys_to_delete:
            del peak_dic[key]

        keys_to_delete = []

        for key in decon_peak_dic:
            if decon_peak_dic[key].empty:
                keys_to_delete.append(key)

        for key in keys_to_delete:
            del decon_peak_dic[key]

        for ion in peak_dic.keys():
            data_df = smooth_df_final[ion].to_frame()
            data_df = data_df.dropna(axis=0, how="any")
            #################
            # 重要：使用to_frame()可将series转换为单列（单行）的df，方便格式统一
            #################
            peak_dic[ion] = self.sharpness(peak_dic[ion], data_df)

        for ion in decon_peak_dic.keys():
            new_peak_df = pd.DataFrame(columns=["left", "apex", "right", "SV"])
            for peak in decon_peak_dic[ion].index.values:
                data_df = decon_data_df.loc[:, [str(peak)]]
                data_df = data_df.dropna(axis=0, how="any")
                # print("debug_csv_data_df = ", data_df)
                peak_df = decon_peak_dic[ion].loc[[peak], :]
                #################
                # 重要：切片单行（单列）时，将目标行（列）外面再加一个[]，可实现切片后保持df的格式
                #################
                peak_df = self.sharpness(peak_df, data_df)
                # print("debug_peak_df = ", peak_df)
                new_peak_df = pd.concat((new_peak_df, peak_df))
                #################
                # 重要：列名一致时，使用df = pd.concat((df1, df2))可纵向拼接两个df，注意将拼接后的df赋值给一个变量
                #################
            # print("debug_new_peak_df = ", new_peak_df)

            decon_peak_dic[ion] = new_peak_df
            # decon_df列名改成解卷积后的apex

        return peak_dic, decon_peak_dic

    def group_peak(self, peak_dic, decon_peak_dic, smooth_df, smooth_df_final, decon_data_df, wid, sigma, bin_num,
                   group_method):
        quant_result_df = pd.DataFrame(columns=['Peak_group', 'Ion', "Relative_Peak_Area", "Peak_Height"])
        # 注意数据中有3套RT，使用时不要混淆：smooth_df raw RT, smooth_df_final 3点校正RT, peak_group_df, 切分RT
        rt_list = smooth_df.index.values.tolist()
        # 将原df的rt的每个时间点之间等分成10份，并加头加尾，形成新rt list
        if bin_num >= 1:
            bin_num = round(bin_num)
            new_list = [rt_list[0]] + [rt_list[i] + (rt_list[i + 1] - rt_list[i]) / bin_num * j for i in
                                       range(len(rt_list) - 1)
                                       for j in range(1, bin_num + 1)] + [rt_list[-1]]
        elif 0 < bin_num < 1:
            bin_num_re = math.ceil(1 / bin_num)
            new_list = [rt_list[0]]  # 新列表中的第一个元素和原列表相同
            idx = bin_num_re  # 下一个需要取的索引
            while idx < len(rt_list) - 1:  # 循环直到快要取到最后一个点
                new_list.append(rt_list[idx])  # 将需要取的元素加入新列表
                idx += bin_num_re  # 计算下一个需要取的索引
            new_list.append(rt_list[-1])  # 新列表中的最后一个元素和原列表相同
        else:
            print("bin_num_error")

        matched_wave = pd.DataFrame(index=new_list, columns=["rt", "SV", "left", "right", "ions"])
        matched_wave["rt"] = matched_wave.index
        matched_wave = matched_wave.fillna(0)
        matched_wave['ions'] = [pd.DataFrame() for i in range(len(matched_wave))]
        for ion in peak_dic.keys():
            for peak, row in peak_dic[ion].iterrows():
                apex = row[1]
                left = row[0]
                right = row[2]
                for i in range(0, (len(matched_wave) - 1)):

                    if matched_wave.index[i] <= apex < matched_wave.index[i + 1]:

                        matched_wave.at[matched_wave.index[i], "SV"] += peak_dic[ion].loc[peak, "SV"]

                        # 找最小的left
                        if matched_wave.iloc[i, 2] == 0:
                            matched_wave.iloc[i, 2] = left
                        elif matched_wave.iloc[i, 2] > left:
                            matched_wave.iloc[i, 2] = left
                        # 找最大的right
                        if matched_wave.iloc[i, 3] == 0:
                            matched_wave.iloc[i, 3] = right
                        elif matched_wave.iloc[i, 3] < right:
                            matched_wave.iloc[i, 3] = right
                            # ion list append当前ion
                        matched_wave.iloc[i, 4].loc[ion, "intensity"] = smooth_df_final.loc[apex, ion]

                        break

        for ion in decon_peak_dic.keys():
            for peak_1, row_1 in decon_peak_dic[ion].iterrows():
                apex_1 = row_1[1]
                left_1 = row_1[0]
                right_1 = row_1[2]
                for i in range(len(matched_wave.index) - 1):
                    if matched_wave.index[i] <= apex_1 < matched_wave.index[i + 1]:

                        matched_wave.at[matched_wave.index[i], "SV"] += decon_peak_dic[ion].loc[peak_1, "SV"]

                        if matched_wave.iloc[i, 2] == 0:
                            matched_wave.iloc[i, 2] = left_1
                        elif matched_wave.iloc[i, 2] > left_1:
                            matched_wave.iloc[i, 2] = left_1

                        if matched_wave.iloc[i, 3] == 0:
                            matched_wave.iloc[i, 3] = right_1
                        elif matched_wave.iloc[i, 3] < right_1:
                            matched_wave.iloc[i, 3] = right_1

                        matched_wave.iloc[i, 4].loc[ion, "intensity"] = decon_data_df.loc[apex_1, peak_1]

                        break
        # print(matched_wave)

        # sg_list = []
        # for n in range(0, wid * 2 + 1):
        #     sg_list.append((1 - (pow(((-wid + n) / sigma), 2))) * (math.exp(-0.5 * pow(((-wid + n) / sigma), 2))))

        # # lambda x:后的函数确认无误
        # matched_wave["gauss_SV"] = matched_wave["SV"].rolling(window=wid * 2 + 1, min_periods=wid * 2 + 1,
        #                                                       center=True).apply(
        #     lambda x: sum(np.multiply(np.asarray(x), np.asarray(sg_list)))).to_list()
        # matched_wave = matched_wave.dropna(axis=0, how="any")
        # matched_wave.plot(kind="line", y=["SV", "gauss_SV"])
        # plt.savefig('./results/group/{}.jpg'.format("matched_wave"))
        # plt.show()

        # # peak_group_df = pd.DataFrame(index=matched_wave.index)
        # peak_group_df = pd.DataFrame(columns=matched_wave.columns)
        # # print("debug_peak_group_df = ", peak_group_df)
        # for n in range(1, matched_wave.shape[0] - 1):
        #     if matched_wave.iloc[n, 5] > 0 and matched_wave.iloc[n, 5] - matched_wave.iloc[n - 1, 5] > 0 and \
        #             matched_wave.iloc[n, 5] - matched_wave.iloc[n + 1, 5] > 0:
        #         peak_group_df.loc[len(peak_group_df)] = matched_wave.iloc[n]

        if group_method == 0:
            matched_wave, peak_group_df = self.MSDIAL(matched_wave, wid, sigma)
        elif group_method == 1:
            matched_wave, peak_group_df = self.AMDIS(matched_wave)

        ##########
        for group, row in peak_group_df.iterrows():
            quant_df = row["ions"]
            quant_df_index = quant_df.index.tolist()
            Quant_Ion = int(quant_df_index[0])
            print(Quant_Ion)
            if row["rt"] in decon_data_df.index.values and decon_data_df.loc[
                row["rt"], lambda d: d.columns.str.contains(f'{Quant_Ion}_peak')].isnull().values.all() == False:
                Quant_Ion_try_list = []
                for Quant_Ion_try in quant_df_index:
                    if (row["rt"] in decon_data_df.index.values and decon_data_df.loc[
                        row["rt"], lambda d: d.columns.str.contains(
                            f'{Quant_Ion_try}_peak')].isnull().values.any()) or (
                            row["rt"] not in decon_data_df.index.values):
                        Quant_Ion_try_list.append(Quant_Ion_try)
                        sliced_df = smooth_df_final.loc[
                            (smooth_df_final.index >= row["left"]) & (smooth_df_final.index <= row["right"]), [
                                Quant_Ion_try]]
                        peak_group_df.loc[group, "Relative_Peak_Area"] = sliced_df.iloc[:, 0].sum()
                        peak_group_df.loc[group, "Peak_Height"] = sliced_df.iloc[:, 0].max()
                        peak_group_df.loc[group, "Quant_Ion"] = int(Quant_Ion_try)
                if not Quant_Ion_try_list:
                    col = int(str(row["rt"]).split('.')[0])
                    print(decon_data_df.filter(like=f"{Quant_Ion}_peak_{col}"))
                    sliced_df = decon_data_df.loc[
                        (decon_data_df.index >= row["left"]) & (decon_data_df.index <= row["right"]),
                            :]
                    sliced_df = sliced_df.filter(like=f"{Quant_Ion}_peak_{col}")  # todo: rt时间不匹配，无法正确识别。
                    if sliced_df.empty == False:
                        peak_group_df.loc[group, "Relative_Peak_Area"] = sliced_df.iloc[:, 0].sum()
                        peak_group_df.loc[group, "Peak_Height"] = sliced_df.iloc[:, 0].max()
                        peak_group_df.loc[group, "Quant_Ion"] = Quant_Ion
                    else:
                        sliced_df = smooth_df_final.loc[
                            (smooth_df_final.index >= row["left"]) & (smooth_df_final.index <= row["right"]), [
                                Quant_Ion]]
                        peak_group_df.loc[group, "Relative_Peak_Area"] = sliced_df.iloc[:, 0].sum()
                        peak_group_df.loc[group, "Peak_Height"] = sliced_df.iloc[:, 0].max()
                        peak_group_df.loc[group, "Quant_Ion"] = int(Quant_Ion)

            else:
                sliced_df = smooth_df_final.loc[
                    (smooth_df_final.index >= row["left"]) & (smooth_df_final.index <= row["right"]), [
                        Quant_Ion]]
                peak_group_df.loc[group, "Relative_Peak_Area"] = sliced_df.iloc[:, 0].sum()
                peak_group_df.loc[group, "Peak_Height"] = sliced_df.iloc[:, 0].max()
                peak_group_df.loc[group, "Quant_Ion"] = int(Quant_Ion)

            for select_ion in quant_df_index:
                quant_result_df_len = len(quant_result_df)
                sliced_df = smooth_df_final.loc[
                    (smooth_df_final.index >= row["left"]) & (smooth_df_final.index <= row["right"]), [
                        select_ion]]
                quant_result_df.loc[quant_result_df_len, "Relative_Peak_Area"] = sliced_df.iloc[:, 0].sum()
                quant_result_df.loc[quant_result_df_len, "Peak_Height"] = sliced_df.iloc[:, 0].max()
                quant_result_df.loc[quant_result_df_len, "Ion"] = int(select_ion)
                quant_result_df.loc[quant_result_df_len, "Peak_group"] = row["rt"]

            # left = row[2]
            # right = row[3]
            # x = smooth_df_final.loc[left: right, :].index.values.tolist()
            # n = 1
            # legend_list = []
            # colormap = plt.cm.gist_ncar

            # for i in peak_group_df.loc[group, "ions"].index.values.tolist():
            #     y = smooth_df_final.loc[left: right, i]
            #     # Ignore missing y-values and plot only the available values
            #     x_valid = [x_val for x_val, y_val in zip(x, y) if not pd.isnull(y_val)]
            #     y_valid = [y_val for y_val in y if not pd.isnull(y_val)]
            #     #line, = plt.plot(x_valid, y_valid, color=colormap(n))
            #     #legend_list.append((line, str(i), peak_group_df.loc[group, "ions"].loc[i, 'intensity']))
            #     # Append a tuple with line, curve label, and maximum value
            #     n = n + 8
            #
            #     # Sort the legend list based on the maximum value of each curve
            # legend_list.sort(key=operator.itemgetter(2), reverse=True)
            #
            # if len(legend_list) > 30:
            #     legend_list = legend_list[:10]  # Keep only the first ten items
            #     legend_lines = [item[0] for item in legend_list]  # Extract the lines from the legend list
            #     legend_labels = [item[1] for item in legend_list]  # Extract the curve labels from the legend list
            #     plt.legend(legend_lines, legend_labels, ncol=math.ceil((len(legend_list)) / 30),
            #                bbox_to_anchor=(1.05, 1),
            #                loc='upper left', borderaxespad=0.)
            # else:
            #     legend_list = legend_list[:10]
            #     legend_lines = [item[0] for item in legend_list]  # Extract the lines from the legend list
            #     legend_labels = [item[1] for item in legend_list]  # Extract the curve labels from the legend list
            #     plt.legend(legend_lines, legend_labels, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
            #
            # plt.gcf().set_size_inches(16, 8)
            # #plt.savefig('./results/group/{}.jpg'.format(str(peak_group_df.loc[group, "rt"])))
            # plt.close()
        # matched_wave.to_excel('./results/group/matched_wave.xlsx')
        # peak_group_df.to_excel('./peak_group_df.xlsx')

        return matched_wave, peak_group_df, quant_result_df

    def peak_group_df_add_retention_infor(self, peak_group_df, RI_Presentation, standard_df, RI_min, RI_max):
        """
        peak_group_df: 峰分组信息
        RI_Presentation: RI的计算模式
        standard_df: RI计算参考标准表
        RI_min, RI_max: RI阈值
        """
        # todo: new def --jfyd
        for i, row in peak_group_df.iterrows():
            rt_sample = row['rt'] / 60
            if RI_Presentation == "Kovats":
                peak_group_df.loc[i, "Kovats_RI"] = self.RT_to_Kovats_RI_transform(rt_sample, standard_df, RI_min,
                                                                                   RI_max)
            elif RI_Presentation == "Fiehn":
                peak_group_df.loc[i, "Fiehn_RI"] = self.RT_to_Fiehn_RI_transform(rt_sample, standard_df)
        return peak_group_df

    def MSDIAL(self, matched_wave, wid, sigma):
        """
        MSDIAL找物质
        """
        sg_list = []
        for n in range(0, wid * 2 + 1):
            sg_list.append((1 - (pow(((-wid + n) / sigma), 2))) * (math.exp(-0.5 * pow(((-wid + n) / sigma), 2))))

        # lambda x:后的函数确认无误
        matched_wave["gauss_SV"] = matched_wave["SV"].rolling(window=wid * 2 + 1, min_periods=wid * 2 + 1,
                                                              center=True).apply(
            lambda x: sum(np.multiply(np.asarray(x), np.asarray(sg_list)))).to_list()
        matched_wave = matched_wave.dropna(axis=0, how="any")
        # matched_wave.plot(kind="line", y=["SV", "gauss_SV"])
        # plt.savefig('./results/group/{}.jpg'.format("matched_wave"))
        # plt.close()

        peak_group_df = pd.DataFrame(columns=matched_wave.columns)
        # print("debug_peak_group_df = ", peak_group_df)
        for n in range(1, matched_wave.shape[0] - 1):
            if matched_wave.iloc[n, 5] > 0 and matched_wave.iloc[n, 5] - matched_wave.iloc[n - 1, 5] > 0 and \
                    matched_wave.iloc[n, 5] - matched_wave.iloc[n + 1, 5] > 0:
                peak_group_df.loc[len(peak_group_df)] = matched_wave.iloc[n]
        peak_group_df = peak_group_df.drop(peak_group_df[peak_group_df['SV'] == 0].index)


        return matched_wave, peak_group_df

    def AMDIS(self, matched_wave):
        """
        AMDIS 算法复现
        """
        matched_wave.insert(5, "compound_point", np.nan, allow_duplicates=False)
        # matched_wave.insert(2, "left_bound", np.nan, allow_duplicates=False)
        # matched_wave.insert(3, "right_bound", np.nan, allow_duplicates=False)

        for i, rows in matched_wave.iterrows():
            sv = rows["SV"]
            if sv != 0:
                index_list = list(matched_wave.index)
                ii = index_list.index(i)
                bound = max(int(500 / matched_wave.iloc[ii - 1:ii + 2, 1].sum()), 1)
                ii_b = max(ii - bound, 0)
                ii_a = min(ii + bound, len(index_list))
                sv_series = matched_wave.iloc[ii_b:ii_a + 1, 1]

                if not sv_series.empty:
                    if sv == max(list(sv_series.values)):
                        matched_wave.loc[i, ["compound_point"]] = sv
                        matched_wave.loc[i, ["left"]] = min(
                            [i for i in matched_wave.iloc[ii_b:ii_a + 1, 2].values if i != 0])
                        matched_wave.loc[i, ["right"]] = max(
                            [i for i in matched_wave.iloc[ii_b:ii_a + 1, 3].values if i != 0])
                        for j in range(ii_b, ii_a + 1):
                            matched_wave.loc[i, ["ions"]][0] = matched_wave.loc[i, ["ions"]][0].combine_first(
                                matched_wave.iloc[j, [4]][0])

        peak_group_df = matched_wave.dropna(axis=0, how="any")
        # matched_wave.plot(kind="line", y=["SV", "compound_point"])
        # plt.savefig('./results/group/{}.jpg'.format("matched_wave"))
        # plt.close()

        return matched_wave, peak_group_df

    ##################

    ##################
    def is_number(self, s):
        """
        判断字符串 s 是否为数字
        """
        pattern = r'^[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?$'
        return bool(re.match(pattern, s))

    def remove_content_from_list(self, lst, content):
        pattern = re.escape(content)  # 转义特殊字符
        regex = re.compile(pattern)
        return [regex.sub('', string) for string in lst if isinstance(string, str)]

    def read_msp(self, msp):
        msp_file = open(msp, "r")
        list_1 = msp_file.readlines()
        # print(list_1)
        new_list = [item.replace('NAME:', 'Name:') for item in list_1]
        list_1 = [item.replace("Num peaks:", "Num Peaks:") for item in new_list]
        # list_1 = replace_name(list_1)
        # print(list_1)
        lists = str(list_1)
        lines = lists.split("Name: ")
        meta = {}
        for l in lines:
            line1 = l.strip().split("\\n")
            # print(line1)
            name_1 = line1[0]
            # print("name = ", name_1)
            line2 = l.strip().split("Num Peaks:")
            ion_intens_dic = {}
            # print("line2 = ", line2)
            if len(line2) > 1:
                if ';' in line2[1]:
                    # 使用正则表达式匹配数字对
                    matches = re.findall(r"(\d+) (\d+);", line2[1])
                    # 创建字典
                    ion_intens_dic = {}
                    for key, value in matches:
                        key = round(float(key))
                        value = int(value)
                        if key in ion_intens_dic:
                            ion_intens_dic[key] = max(ion_intens_dic[key], value)
                        else:
                            ion_intens_dic[key] = value

                elif '\\t' in line2[1]:
                    line3 = line2[1].split("\\n', '")[1:-2]
                    for ion in line3:
                        ion1 = ion.split("\\t")
                        # print("ion1 = ", ion1)
                        if len(ion1) == 2 and self.is_number(ion1[0]) and self.is_number(ion1[1]):
                            key = round(float(ion1[0]))
                            value = float(ion1[1])
                            if key in ion_intens_dic:
                                ion_intens_dic[key] = max(ion_intens_dic[key], value)
                            else:
                                ion_intens_dic[key] = value
                elif '\\n' in line2[1]:
                    line3 = line2[1].split("\\n', '")[1:-2]
                    for ion in line3:
                        ion1 = ion.split(" ")
                        # print("ion1 = ", ion1)
                        if len(ion1) == 2 and self.is_number(ion1[0]) and self.is_number(ion1[1]):
                            key = round(float(ion1[0]))
                            value = float(ion1[1])
                            if key in ion_intens_dic:
                                ion_intens_dic[key] = max(ion_intens_dic[key], value)
                            else:
                                ion_intens_dic[key] = value
                else:
                    print('格式无法识别')
                    # print("ion_intens_dic = ", ion_intens_dic)
            meta[name_1] = ion_intens_dic
        RI_df = self.read_msp_RI(msp)
        return meta, RI_df

    def group_cmp_inf(self, lines):
        """
        FUNCTION:检索物质名所在行
        """
        group_inf_idx = [i for i, p in enumerate(lines) if 'Name:' in p]
        group_inf_idx.append(len(lines))
        return group_inf_idx

    def read_msp_RI(self, msp):
        """
        FUNCTION:转换样本RT-RI
        RT_data: RT库
        msp: MSP库
        return: msp_RI
        """
        msp_file = open(msp, "r")
        lines = msp_file.readlines()
        group_inf_idx = self.group_cmp_inf(lines)
        RI_df = pd.DataFrame(columns=['Name', 'RI'])
        for j in range(len(group_inf_idx) - 1):
            group_inf = lines[group_inf_idx[j]:group_inf_idx[j + 1]]
            # 定义要匹配的字符串前缀
            prefixes = [r'SemiStdNP=\d+', r'RI:\d+\n']
            # 定义正则表达式
            pattern = "|".join(prefixes)
            for string in group_inf:
                if 'Name:' in string:
                    RI_list = [string.replace('Name: ', '')]
                if matches := re.findall(pattern, string):
                    RI_list.extend([int(re.findall(r"\d+", match)[0]) for match in matches])
                    # print(RI_list)
                    RI_df.loc[len(RI_df.index)] = RI_list
        RI_df['Name'] = RI_df['Name'].str.rstrip('\n')
        RI_df.set_index('Name', inplace=True)
        return RI_df

    ########
    # backup0215中第1部分原位置
    ########

    def dot_product_distance(self, p, q):
        if (np.sum(p)) == 0 or (np.sum(q)) == 0:
            score = 0
        else:
            score = np.power(np.sum(q * p), 2) / \
                    (np.sum(np.power(q, 2)) * np.sum(np.power(p, 2)))
        # return 1 - np.sqrt(score)
        # 由于写的是Dot product distance，即距离，因此先直接导出score，与NIST原文一致
        return score

    def weighted_dot_product_distance(self, compare_df):
        # print("debug compare_df = ")
        # print(compare_df)
        m_q = pd.Series(compare_df.index)
        # m_q是df的index，是输入的离子
        # print("debug m_q = ")
        # print(m_q)
        i_q = np.array(compare_df.iloc[:, 0])
        # i_q是df的第一列，是实测响应
        # 后续i_q与m_q要相乘，只有array或Series才可以相乘，因此要把df转为array
        # print("debug i_q = ")
        # print(i_q)
        i_r = np.array(compare_df.iloc[:, 1])
        # i_r是df的第二列，是lib响应
        # print("debug i_r = ")
        # print(i_r)
        k = 0.5
        # NIST的k=0.5，改为1可提高丰度比的权重
        l = 2
        w_q = np.power(i_q, k) * np.power(m_q, l)
        # print("debug w_q = ")
        # print(w_q)
        # print(type(w_q))
        w_r = np.power(i_r, k) * np.power(m_q, l)
        # print("debug w_r = ")
        # print(w_r)

        # 如果组某离子或所有离子读数为0，该步也可正常计算
        ss = self.dot_product_distance(w_q, w_r)
        # print("debug ss = ", ss)
        shared_spec = np.vstack((i_q, i_r))
        shared_spec = pd.DataFrame(shared_spec)
        shared_spec = shared_spec.loc[:, (shared_spec != 0).all(axis=0)]
        # print("debug_shared_spec = ", shared_spec)
        # 取共有离子
        m = int(shared_spec.shape[1])
        # print("debug_m = ", m)
        # 如果要提高丰度比的权重，即增加m，那么要在该处：composite_score = ((NU*ss) + (m*ave_FR)) / (NU + m)增加m，因为下面有个m是否大于等于2的判断
        if m >= 2:
            FR = 0
            for i in range(1, m):
                # df.iat中，行数在前，列数在后，取值时0是第一行/列，因此range是1到n
                # range中1, n,包含1，但不包含n，最大n-1
                s = (shared_spec.iat[0, i] / shared_spec.iat[0, (i - 1)]) * (
                        shared_spec.iat[1, (i - 1)] / shared_spec.iat[1, i])
                if s > 1:
                    s = 1 / s
                FR = FR + s
            ave_FR = FR / (m - 1)
            NU = int(len(compare_df))
            # 原array中行数是物质包含的离子数
            composite_score = ((NU * ss) + (m * ave_FR)) / (NU + m)
        else:
            composite_score = ss

        return composite_score

    ########
    # backup0215中第2部分原位置
    ########

    ############

    def euclidean_intensity(self, x, y, k, b):
        '''y = kx + b '''
        return math.fabs(b + k * (x - 0) - y) / (math.sqrt(k * k + 1))

    def line_between_dots(self, x1, y1, x2, y2):
        '''y = kx + b 返回 k,b'''
        if x1 == x2:
            return f'input_error: x1 = x2 = {x1}', np.nan
        else:
            return (y1 - y2) / (x1 - x2), (x1 * y2 - x2 * y1) / (x1 - x2)

    def line(self, x, k, b):  # 定义拟合函数形式
        return k * x + b

    def base_line(self, rt_list, intensity_list):
        para, pcov = curve_fit(self.line, rt_list, intensity_list)
        y_fitted = self.line(rt_list, para[0], para[1])  # 画出拟合后的曲线
        # plt.plot(rt_list, y_fitted, '-b')
        return para[0], para[1], y_fitted

    def baseline_correction(self, peak_df, smooth_df, ion, smooth_df_final):
        base_line_rt = []
        base_line_intensity = []
        list_now = []
        for peak_line in peak_df.itertuples():
            # print(peak_line)
            left = getattr(peak_line, 'left')
            right = getattr(peak_line, 'right')
            apex = getattr(peak_line, 'apex')
            son_df = smooth_df.loc[left:right, [ion]]
            son_df_l = smooth_df.loc[left:apex, [ion]]
            son_df_r = smooth_df.loc[apex:right, [ion]]
            min_l_dots = son_df_l.nsmallest(1, ion, keep='all')
            min_r_dots = son_df_r.nsmallest(1, ion, keep='all')
            min2dots = pd.concat([min_l_dots, min_r_dots], axis=0)
            k, b = self.line_between_dots(min2dots.index[0], min2dots.iloc[0, 0], min2dots.index[1],
                                          min2dots.iloc[1, 0])
            euclidean_intensity_list = []
            for rt, intensity in son_df.iterrows():
                euclidean_intensity_list.append(intensity.values[0] - self.line(rt, k, b))
                # 此处计算的是原丰度减临时基线y值得到的丰度，而不计算点到线的垂直距离，因为根据前者计算完后，顶点仍是顶点，且峰型不会有太大变化
            son_df["euclidean_intensity_list"] = euclidean_intensity_list
            list_now.extend(list(son_df[ion]))
            son_son_df = son_df.nsmallest(int(son_df.shape[0] / 2), "euclidean_intensity_list", keep='all')  # 取了最小的一半df
            ##计算全部点到基线的距离
            k1, b1, y_fitted = self.base_line(son_son_df.index, son_son_df.iloc[:, 0])
            for rt1, intensity1 in son_df.iterrows():
                new_intensity = intensity1.values[0] - self.line(rt1, k1, b1)
                if new_intensity >= 0:
                    base_line_intensity.append(new_intensity)
                else:
                    base_line_intensity.append(0)
                base_line_rt.append(rt1)
            # plt.plot(base_line_rt, list_now, 'o')
            # plt.ylim(min(son_df[ion]) * 0.8, max(son_df[ion]) * 1.2)
            # plt.xlim(min(son_df.index) - 5, max(son_df.index) + 5)
            # plt.savefig('./results/peak/{}.png'.format(str(ion) + "_" + str(apex) + "_adjusted_baseline"))
            # plt.show()
        for index1, line1 in smooth_df.iterrows():
            if index1 in base_line_rt:
                smooth_df_final.loc[index1, ion] = base_line_intensity[base_line_rt.index(index1)]

        return smooth_df_final

    ############
    def pre_de_redundancy_peak(self, smooth_df, ion_list):
        """
        去除冗余峰
        1. 先从smooth_df中取出[ion]dataframe，再循环13个取
        2. if anyone = 0 -> ?
        3. if anyone != 0 -> number of >mean dots > 7 -> median_dot
        """
        all_NF_list = []
        for ion in ion_list:
            m = 0
            i_list = smooth_df[ion].values.tolist()
            NF_list = []
            while m < len(i_list):
                window_list = i_list[m:m + 13]
                if 0 not in window_list:
                    window_mean = np.mean(window_list)
                    number = sum(i > window_mean for i in window_list)
                    # print("debug_number = ", number)
                    if number >= 7:
                        median_value = median(window_list)
                        NF_list.append(abs(median_value - window_mean) / (window_mean ** 0.5))
                m = m + 13
                # print("debug_m = ", m)
            all_NF_list.extend(NF_list)
        NF = median(all_NF_list)
        print('Noise_factor = ', NF)

        return NF

    def de_redundancy_peak(self, smooth_df_final, smooth_df, peak_df, NF, ion):
        index_list = []
        for index, row in peak_df.iterrows():
            # print("debug index = ", index)
            tmp_apex_intensity = smooth_df.loc[row["apex"], ion]
            ex_tmp_apex_intensity = smooth_df_final.loc[row["apex"], ion]
            if (4 * NF * ((tmp_apex_intensity) ** 0.5)) < ex_tmp_apex_intensity:
                index_list.append(index)
        peak_df.drop(index_list, axis=0)

        return peak_df

    ############################

    ############################
    def RT_similarity(self, rt_s, rt_l, window, level_factor, max_penalty):
        """
        NIST中采用的比较RI相似性方式
        rt_s:样品rt.
        rt_l：库中rt.
        penalty_rate：罚分程度. weak-20,average-50,strong-100,very strong-200
        rt_tolerance：rt容忍度.
        """
        delta_rt = abs(rt_s - rt_l)
        if delta_rt < window:
            rt_penalty = 0
        else:
            rt_penalty = ((delta_rt / window) - 1) * level_factor
            if rt_penalty > max_penalty:
                rt_penalty = max_penalty

        return rt_penalty, delta_rt

    def RI_similarity(self, ri_s, ri_l, ri_window, ri_window_scale, level_factor, max_penalty):
        # default: ri_s, ri_l, 40, 5, 0.05, 0.1

        delta_ri = abs(ri_s - ri_l)

        if ri_window_scale > 0:
            ri_window = ri_window + ri_window_scale * 0.01 * ri_s

        if delta_ri <= ri_window:
            ri_penalty = 0
        else:
            ri_penalty = ((delta_ri / ri_window) - 1) * level_factor
            if ri_penalty > max_penalty:
                ri_penalty = max_penalty

        return ri_penalty, delta_ri

    def New_Compound(self, group_rt_s, row, score_df, results_df, min_ions_number: int, limit_of_score: float):
        # yuannote: New_Compound函数中，min_ions_number对应不使用-向导页面3、RT-向导页面3、RI-向导页面3中，新物质识别参数中
        # 参与新物质识别最少离子数 对应的数字
        # limit_of_score对应上述3个页面中，新物质相似性阈值 对应的数字

        """
        提示可能为新化合物
        row: peak_group_df 遍历的行信息
        score_df：库中化合物匹配得分df
        min_ions_number： 该信号所含最小离子数
        return: 注释了新化合物的results_df
        """
        if len(row["ions"]) >= min_ions_number:
            if max(score_df["score"]) < limit_of_score:
                results_df.loc[group_rt_s, "New_Compound"] = "这个信号对应的化合物可能暂时没有包含在提供的参考库中"
        return results_df

    def species_pool_opt(self, results_df, species_pool, species_pool_coefficient=0.3,
                         all_match_optimization_coefficient=0.2):
        """
        根据输入的物种库信息优化定性结果
        """
        for result in results_df.itertuples():
            range_left = result[0] / 60 - species_pool_coefficient
            range_right = result[0] / 60 + species_pool_coefficient
            if not result[2].empty:
                optimization_df = species_pool[(species_pool['RT'] < range_right) & (species_pool['RT'] > range_left)]
                all_match_optimization_df = result[2][
                    (result[2]['Score'] > (max(result[2]['Score']) - all_match_optimization_coefficient))]
                IntersectionSet = list(
                    set(all_match_optimization_df.index.tolist()) & set(optimization_df['Name'].tolist()))
                if len(IntersectionSet) == 1:
                    results_df.loc[result[0], "species_pool_alignment_optimization"] = IntersectionSet[0]

        return results_df

    def mass_spectra_comparison_plot(self, group_rt_s, compound, compare_df_1):
        """
        FUNCTION：生成质谱对比图
        group_rt_s: 分组名称
        compound:库中化合物名称
        compare_df_1：包含两个质谱图信息的dataframe
        """

        compare_df_new = compare_df_1[compare_df_1 >= 10]
        compare_df_new.dropna(inplace=True)

        # if not compare_df_new.empty:
        #     # 设置每个点的 x 坐标
        #     x = compare_df_new.index
        #
        #     # 计算图形尺寸和横坐标间距
        #     num_data_points = max(x) / 6
        #     fig_width = max(8, num_data_points)  # 根据数据点数量设置图形宽度
        #
        #     # 创建图形对象
        #     fig, ax = plt.subplots(figsize=(fig_width, 6), dpi=300)
        #
        #     # 去掉外边框
        #     ax.axis('off')
        #
        #     # 绘制垂直线状棒状图
        #     ax.vlines(x, ymin=0, ymax=compare_df_new['intensity_s'], color='r')
        #     ax.vlines(x, ymin=0, ymax=-compare_df_new['intensity_l'], color='b')
        #
        #     ax.set_xlabel('Mass')
        #     ax.set_ylabel('Intensity')
        #
        #     # 在第一个质谱数据集的数据点旁边添加横坐标标注
        #     for i, (mass, intensity) in enumerate(zip(compare_df_new.index, compare_df_new['intensity_s'])):
        #         if intensity > 20:
        #             ax.annotate(round(mass, 2), xy=(mass, intensity), xytext=(0, 5), textcoords='offset points',
        #                         ha='center', va='bottom', fontsize=8)
        #
        #     # 在第二个质谱数据集的数据点旁边添加横坐标标注
        #     for i, (mass, intensity) in enumerate(zip(compare_df_new.index, -compare_df_new['intensity_l'])):
        #         if intensity < -20:
        #             ax.annotate(round(mass, 2), xy=(mass, intensity), xytext=(0, -5), textcoords='offset points',
        #                         ha='center', va='top', fontsize=8)
        #
        #     # 绘制水平线
        #     ax.axhline(y=0, color='black', linewidth=1)
        #
        #     # 设置坐标轴范围
        #     ax.set_ylim(-999, 999)
        #     x_min = min(x)
        #     x_max = max(x)
        #     ax.set_xlim(x_min - 2, x_max + 2)
        #
        #     # 自动调整标签的间距
        #     plt.tight_layout()
        #
        #     compound_name = self.validateTitle(str(compound))
        #
        #     #plt.savefig('./results/comparison_plot/{}.png'.format(str(group_rt_s) + "_" + compound_name))
        #     plt.close()

    def peak_group_identification(self, peak_group_df,
                                  smooth_df,
                                  smooth_df_final,
                                  meta,
                                  RT_df,
                                  decon_peak_dic,
                                  decon_data_df,
                                  peak_group_search,
                                  match_weight,
                                  r_match_weight,
                                  search_wid,
                                  retention_score_mode,
                                  window,
                                  window_scale,
                                  level_factor,
                                  max_penalty,
                                  no_info_penalty,
                                  minimum_number_of_ions,
                                  sim_threshold,
                                  group_weight,
                                  direct_weight,
                                  inaccurate_ri_threshold,
                                  inaccurate_ri_level_factor
                                  ):

        results_df = pd.DataFrame(index=peak_group_df["rt"].values.tolist(),
                                  columns=["Best_match_name", "Best_match", "All_match", "All_match_list", "Quant_Ion",
                                           "Relative_Peak_Area", "Peak_Height", "New_Compound"])
        results_df["Best_match"] = [pd.DataFrame(columns=["Score"]) for i in range(len(results_df))]
        results_df["All_match"] = [pd.DataFrame(columns=["Score"]) for i in range(len(results_df))]
        results_df["All_match_list"] = [[] for i in range(len(results_df))]

        for group, row in peak_group_df.iterrows():

            score_df = pd.DataFrame()

            results_list = []

            group_rt_min = (float(peak_group_df.loc[group, "rt"])) / 60
            group_rt_s = float(peak_group_df.loc[group, "rt"])
            # print("group_rt_min = ", group_rt_min)
            # print("group_rt_s = ", group_rt_s)
            # print(group_rt_min - search_wid)
            if retention_score_mode == "None":
                # print("right = ", right_idx)
                compound_list = RT_df.index.values.tolist()
            elif retention_score_mode == "RI":
                if 'Kovats_RI' in peak_group_df.columns:
                    if not isinstance(peak_group_df.loc[group, "Kovats_RI"], str):
                        group_ri = float(peak_group_df.loc[group, "Kovats_RI"])
                elif 'Fiehn_RI' in peak_group_df.columns:
                    if not isinstance(peak_group_df.loc[group, "Fiehn_RI"], str):
                        group_ri = float(peak_group_df.loc[group, "Fiehn_RI"])
                left_idx = self.find_nearest(RT_df["RI"].values.tolist(), group_ri - search_wid)
                # print("left = ", left_idx)
                # print(group_rt_min + search_wid)
                right_idx = self.find_nearest(RT_df["RI"].values.tolist(), group_ri + search_wid)
                # print("right = ", right_idx)
                compound_list = RT_df.iloc[left_idx: right_idx, :].index.values.tolist()
            elif retention_score_mode == "RT":
                left_idx = self.find_nearest(RT_df["RT"].values.tolist(), group_rt_min - search_wid)
                # print("left = ", left_idx)
                # print(group_rt_min + search_wid)
                right_idx = self.find_nearest(RT_df["RT"].values.tolist(), group_rt_min + search_wid)
                # print("right = ", right_idx)
                compound_list = RT_df.iloc[left_idx: right_idx, :].index.values.tolist()
            # print("compound_list_start = ", compound_list)
            rt_index = self.find_nearest(smooth_df.index.values.tolist(), group_rt_s)
            group_rt_1 = smooth_df.index.values.tolist()[rt_index]
            # print("compound_list = ", compound_list)
            # 共有3组离子：SIM测的，库里物质有的，group有的
            # shared_ion取库离子中，SIM测到的，再与group比较
            sim_ion = []
            for ion in smooth_df.columns.values.tolist():
                if smooth_df.loc[group_rt_1, ion] != 0:
                    if str(ion).isdigit() == True:
                        sim_ion.append(ion)
            # print(group, "_sim_ion = ", sim_ion)

            if peak_group_search == 1:
                compare_df = row["ions"]
                # else:
                n = self.find_nearest(smooth_df.index.values.tolist(), group_rt_s)
                direct_compare_df = pd.DataFrame(columns=['intensity'])
                selected_cols = []
                selected_values = []
                for col_name, value in smooth_df.iloc[n].items():
                    if str(col_name).isdigit() and value > 0:
                        selected_cols.append(col_name)
                        selected_values.append(value)
                direct_compare_df['intensity'] = selected_values
                direct_compare_df.index = selected_cols

            # 归一化处理
            compare_df["intensity"] = compare_df["intensity"].apply(lambda x: (x / compare_df["intensity"].max()) * 999)
            direct_compare_df["intensity"] = direct_compare_df["intensity"].apply(
                lambda x: (x / direct_compare_df["intensity"].max()) * 999)

            # print("compound_list = ", compound_list)

            if len(compound_list) > 0:
                for compound in compound_list:
                    compound_ion = []
                    if compound in meta:
                        for ion_1, intensity in meta[compound].items():
                            if isinstance(intensity, (int, float, np.float64)) or (
                                    isinstance(intensity, str) and intensity.isdigit()):
                                # print(intensity)
                                if int(intensity) > 0:
                                    compound_ion.append(int(ion_1))
                            else:
                                print(compound, "格式有误")
                    # print(compound, "_compound_ion = ", compound_ion)
                    shared_ion = [i for i in compound_ion if i in sim_ion]
                    # shared_ion_1 = [i for i in shared_ion if i in group_ion]

                    # print(compound, "_shared_ion = ", shared_ion)

                    com_df_l_dic = {}
                    for i in shared_ion:
                        com_df_l_dic[i] = meta[compound][int(i)]
                    com_df_l = pd.Series(com_df_l_dic, dtype=float)
                    # print(compound, "com_df_l = ", com_df_l)

                    compare_df_1 = pd.concat([compare_df, com_df_l], axis=1, ignore_index=False)
                    # compare_df_1 = compare_df_1.dropna(axis=0, how="any")
                    compare_df_1.fillna(0, inplace=True)
                    compare_df_1 = compare_df_1.loc[~(compare_df_1 == 0).all(axis=1)]
                    compare_df_1.rename(columns={'intensity': 'intensity_s', 0: 'intensity_l'},
                                        inplace=True)  # inplace = True，表示在原始dataframe上修改列名
                    # compound_name = validateTitle(compound)
                    # compare_df_1.to_excel('./results/score/{}.xlsx'.format(str(group_rt_s) + "_" + str(compound_name)))
                    # print('compare_df_1', compare_df_1)

                    r_match_compare_df = compare_df_1.loc[compare_df_1.iloc[:, 1] != 0]
                    # r_match_compare_df.to_excel('./results/score/{}.xlsx'.format(str(group_rt_s) + "_" + str(compound_name) + "_rmatch"))

                    direct_compare_df_1 = pd.concat([direct_compare_df, com_df_l], axis=1, ignore_index=False)
                    direct_compare_df_1.fillna(0, inplace=True)
                    direct_compare_df_1.rename(columns={'intensity': 'intensity_s_direct', 0: 'intensity_l_direct'},
                                               inplace=True)  # inplace = True，表示在原始dataframe上修改列名
                    R_direct_compare_df_1 = direct_compare_df_1.loc[direct_compare_df_1.iloc[:, 1] != 0]

                    MS_score = self.weighted_dot_product_distance(compare_df_1)
                    R_MS_score = self.weighted_dot_product_distance(r_match_compare_df)
                    R_direct_MS_score = self.weighted_dot_product_distance(R_direct_compare_df_1)
                    score_df.loc[compound, "MS_score"] = MS_score
                    score_df.loc[compound, "R_MS_score"] = R_MS_score
                    score_df.loc[compound, "R_direct_score"] = R_direct_MS_score

                    if retention_score_mode == "None":
                        Retention_score = 0
                    elif retention_score_mode == "RI":
                        if 'Kovats_RI' or 'Fiehn_RI' in peak_group_df.columns:
                            if 'Kovats_RI' in peak_group_df.columns:
                                ri_s = row["Kovats_RI"]
                                # print("ri_s = ", ri_s)
                            elif 'Fiehn_RI' in peak_group_df.columns:
                                ri_s = row['Fiehn_RI']

                            ri_l = RT_df.loc[compound, 'RI']
                            # print("ri_l = ", ri_l)
                            if type(ri_l) not in [float, int, np.float64, np.int64] or ri_l == np.nan:
                                # print(type(ri_l))
                                Retention_score = no_info_penalty
                            elif ri_s <= inaccurate_ri_threshold:
                                Retention_score, delta_ri = self.RI_similarity(ri_s,
                                                                               ri_l,
                                                                               window,
                                                                               window_scale,
                                                                               inaccurate_ri_level_factor,
                                                                               max_penalty)
                                score_df.loc[compound, "Delta_RI"] = delta_ri
                                score_df.loc[compound, "Reference_RI"] = ri_l
                            else:
                                Retention_score, delta_ri = self.RI_similarity(ri_s,
                                                                               ri_l,
                                                                               window,
                                                                               window_scale,
                                                                               level_factor,
                                                                               max_penalty)
                                score_df.loc[compound, "Delta_RI"] = delta_ri
                                score_df.loc[compound, "Reference_RI"] = ri_l

                        else:
                            # todo:写入错误日志，说明没计算物质RI
                            Retention_score = 0

                    elif retention_score_mode == "RT":
                        rt_s = group_rt_min

                        if 'RT' in RT_df.columns:
                            rt_l = RT_df.loc[compound, 'RT']
                            if type(rt_l) not in [float, int, np.float64, np.int64] or rt_l == np.nan:
                                Retention_score = no_info_penalty
                            else:
                                Retention_score, delta_rt = self.RT_similarity(rt_s, rt_l, window, level_factor,
                                                                               max_penalty)
                                score_df.loc[compound, "Delta_RT"] = delta_rt
                                score_df.loc[compound, "Reference_RT"] = rt_l

                        else:
                            # todo: 写入错误日志，说明没有library RT
                            Retention_score = 0

                    score_df.loc[compound, "Retention_score"] = Retention_score

                    if match_weight + r_match_weight == 0 or group_weight + direct_weight == 0:

                        print("score_weight_error")

                        break

                    elif MS_score * R_MS_score * R_direct_MS_score == 0:
                        # 如果MS_score或R_MS或R_direct_MS_score任意一个为0，那么score直接等于0，这个是根据score_df比对非靶结果得出的经验
                        score = 0

                    elif len(row["ions"]) == 1:

                        score = R_direct_MS_score - Retention_score

                    elif len(row["ions"]) == 2:
                        if MS_score <= 0.2:
                            score = (MS_score * 0.1 + R_direct_MS_score * 0.9) - Retention_score

                        elif R_MS_score <= 0.2:
                            score = (R_MS_score * 0.1 + R_direct_MS_score * 0.9) - Retention_score

                        else:
                            group_score = (match_weight * MS_score + r_match_weight * R_MS_score) / (
                                    match_weight + r_match_weight)
                            score = ((group_score * 0.05 + R_direct_MS_score * 0.95) / (
                                    group_weight + direct_weight)) - Retention_score

                    else:

                        if MS_score <= 0.2:

                            score = (MS_score * 0.9 + R_direct_MS_score * 0.1) - Retention_score

                        elif R_MS_score <= 0.2:

                            score = (R_MS_score * 0.9 + R_direct_MS_score * 0.1) - Retention_score

                            # 如果MS_score或R_MS_score过低，则使其更大程度降低最终打分，这个是根据score_df比对非靶结果得出的经验

                        else:
                            group_score = (match_weight * MS_score + r_match_weight * R_MS_score) / (
                                    match_weight + r_match_weight)
                            score = ((group_score * group_weight + R_direct_MS_score * direct_weight) / (
                                    group_weight + direct_weight)) - Retention_score

                    # Rmatch和match得分权重此处修改

                    if score < 0:
                        score = 0

                    score_df.loc[compound, "score"] = score
                    if len(row["ions"]) >= minimum_number_of_ions and score >= sim_threshold:
                        # mass_spectra_comparison_plot(group_rt_s, compound, compare_df_1)
                        # yuannote: 不需要看质谱图比对时，注释掉以加快运行速度
                        results_df.loc[group_rt_s, "All_match"].loc[compound, "Score"] = score
                        results_df.loc[group_rt_s, "All_match"].loc[compound, "mz"] = str(list(compare_df_1.index))
                        results_df.loc[group_rt_s, "All_match"].loc[compound, "intensity_s"] = str(
                            list(compare_df_1['intensity_s']))
                        results_df.loc[group_rt_s, "All_match"].loc[compound, "intensity_l"] = str(
                            list(compare_df_1['intensity_l']))
                        results_df.loc[group_rt_s, "All_match"].loc[compound, "intensity_s_direct"] = str(
                            list(direct_compare_df_1['intensity_s_direct']))
                        results_df.loc[group_rt_s, "All_match"].loc[compound, "intensity_l_direct"] = str(
                            list(direct_compare_df_1['intensity_l_direct']))
                        if retention_score_mode == "RI":
                            if 'Reference_RI' in score_df.columns:
                                results_df.loc[group_rt_s, "All_match"].at[compound, "Reference_RI"] = score_df.loc[
                                    compound, "Reference_RI"]
                            if 'Delta_RI' in score_df.columns:
                                results_df.loc[group_rt_s, "All_match"].at[compound, "Delta_RI"] = score_df.loc[
                                    compound, "Delta_RI"]
                        elif retention_score_mode == "RT":
                            if 'Reference_RT' in score_df.columns:
                                results_df.loc[group_rt_s, "All_match"].at[compound, "Reference_RT"] = score_df.loc[
                                    compound, "Reference_RT"]
                            if 'Delta_RT' in score_df.columns:
                                results_df.loc[group_rt_s, "All_match"].at[compound, "Delta_RT"] = score_df.loc[
                                    compound, "Delta_RT"]

                        results_list.append(compound)
                        score = round(score, 2)
                        results_list.append(score)

            # list作为单独元素，加到df单元格中的时候，双重嵌套即可
            results_df.loc[group_rt_s, "All_match_list"] = [[results_list]]

            # score_df.to_excel('./results/score/{}.xlsx'.format(str(group_rt_s) + "_score_df"))

            print(group_rt_s, "Done")
            # print(results_list)
            # 取所有结果中得分最高的结果
            if results_df.loc[group_rt_s, "All_match"].shape[0] > 0:
                results_df.loc[group_rt_s, "All_match"].sort_values(by="Score", inplace=True, ascending=False)
                name = results_df.loc[group_rt_s, "All_match"].index.values.tolist()[0]
                results_df.loc[group_rt_s, "Best_match"].loc[name, "Score"] = \
                    results_df.loc[group_rt_s, "All_match"].iloc[
                        0, 0]
            elif results_df.loc[group_rt_s, "All_match"].shape[0] == 0 and len(row["ions"]) > 2:
                results_df.loc[group_rt_s, "Best_match_name"] = "Unknown"
            if len(score_df) > 0:
                results_df = self.New_Compound(group_rt_s, row, score_df, results_df,
                                               minimum_number_of_ions, sim_threshold)
            elif len(score_df) == 0:
                results_df.loc[group_rt_s, "New_Compound"] = "这个信号对应的化合物可能暂时没有包含在提供的参考库中"
                # yuannote: 3是3种向导页面3里，新物质识别参数里，参与新物质识别最少离子数，0.1是3种向导页面3里，新物质识别参数里新物质相似性阈值
            # 加入峰面积
            if len(row["ions"]) >= minimum_number_of_ions:
                results_df.loc[group_rt_s, "Relative_Peak_Area"] = peak_group_df.loc[group,
                "Relative_Peak_Area"]
                results_df.loc[group_rt_s, "Peak_Height"] = peak_group_df.loc[group, "Peak_Height"]
                results_df.loc[group_rt_s, "Quant_Ion"] = peak_group_df.loc[group, "Quant_Ion"]
        print(group_rt_s, "Done")
        # 不同group结果相同时，取得分最高的结果，其余结果记为XXX_analogue
        results_df_filtered = pd.DataFrame(index=results_df.index, columns=["Name", "Score"])
        for index, row1 in results_df.iterrows():
            if row1["Best_match"].shape[0] > 0:
                results_df_filtered.loc[index, "Name"] = row1["Best_match"].index.values[0]
                results_df_filtered.loc[index, "Score"] = row1["Best_match"].iloc[0, 0]
        results_df_filtered.sort_values(by='Score', inplace=True, ascending=False)
        results_df_filtered.dropna(axis=0, how="any", inplace=True)
        index_list_filter = results_df_filtered.index.values.tolist()
        # print(results_df_filtered)

        j_list = results_df_filtered["Name"].unique().tolist()
        for j in j_list:
            index_list = [i for i in index_list_filter if results_df_filtered.loc[i, "Name"] == j]
            index_list = index_list[1:]
            for i in index_list:
                name = results_df_filtered.loc[i, "Name"]
                new_name = name + "_analogue"
                results_df.loc[i, "Best_match"].loc[new_name, "Score"] = results_df_filtered.loc[i, "Score"]
                results_df.loc[i, "Best_match"].drop(name, axis=0, inplace=True)
        for index, row in results_df.iterrows():
            if row["Best_match"].shape[0] > 0:
                row["Best_match_name"] = row["Best_match"].index.values.tolist()[0]

        results_df.dropna(subset=["Best_match_name"], inplace=True)
        results_df.drop(results_df.columns[1], axis=1, inplace=True)
        # if not species_pool.empty:
        #     results_df = species_pool_opt(results_df, species_pool)
        #results_df.to_excel('./final_results.xlsx')

        return results_df

    def RT_to_Kovats_RI_transform(self, rt_sample, standard_df, RI_lower_limit, RI_upper_limit):
        """
        FUNCTION:转换样本RT-RI
        rt_sample：样本RT
        standard_df：标品RT-RI
        return:样本RI
        """
        if RI_lower_limit is None:
            RI_lower_limit = 0
        if RI_upper_limit is None:
            RI_upper_limit = 3000
        if np.isnan(rt_sample):
            return "未检索到实测rt内容"
        prev_rows = standard_df.loc[(standard_df['RT (min)'] < rt_sample)].tail(1)
        next_rows = standard_df.loc[(standard_df['RT (min)'] >= rt_sample)].head(1)
        if prev_rows.shape[0] == 0:
            df_sort = standard_df.sort_values('RT (min)', ascending=True).head(2)
        elif next_rows.shape[0] == 0:
            df_sort = standard_df.sort_values('RT (min)', ascending=True).tail(2)
        else:
            df_sort = pd.concat([prev_rows, next_rows])
        RI_low = min(df_sort['RI'])
        RI_high = max(df_sort['RI'])
        RT_low = min(df_sort['RT (min)'])
        RT_high = max(df_sort['RT (min)'])
        ri_sample = round(RI_low + (RI_high - RI_low) * (rt_sample - RT_low) / (RT_high - RT_low))
        if ri_sample < RI_lower_limit:
            return RI_lower_limit
        elif ri_sample > RI_upper_limit:
            return RI_upper_limit
        else:
            return ri_sample

    def RT_to_Fiehn_RI_transform(self, rt_sample, standard_df, Fiehn_RI_list=None):
        if Fiehn_RI_list is None:
            Fiehn_RI_list = [262320, 323120, 381020, 487220, 582620, 668720, 747420, 819620, 886620, 948820, 1006900,
                             1061990, 1113100]
        C8_RT = standard_df['RT (min)'][0] * 60 * 1000
        C8_RI = Fiehn_RI_list[0]
        C9_RT = standard_df['RT (min)'][1] * 60 * 1000
        C9_RI = Fiehn_RI_list[1]
        C28_RT = standard_df['RT (min)'][11] * 60 * 1000
        C28_RI = Fiehn_RI_list[11]
        C30_RT = standard_df['RT (min)'][12] * 60 * 1000
        C30_RI = Fiehn_RI_list[12]
        if rt_sample <= standard_df['RT (min)'][1]:
            z1 = self.linear_fit(C8_RT, C9_RT, C8_RI, C9_RI)
        elif rt_sample >= standard_df['RT (min)'][11]:
            z1 = self.linear_fit(C28_RT, C30_RT, C28_RI, C30_RI)
        else:
            x = np.array(standard_df['RT (min)'][1:12] * 60 * 1000)
            y = np.array(Fiehn_RI_list[1:12])
            z1 = np.polyfit(x, y, 5)  # 用5次多项式拟合
        return np.polyval(z1, rt_sample * 60 * 1000)

    def linear_fit(self, arg0, arg1, arg2, arg3):
        x = np.array([arg0, arg1])
        y = np.array([arg2, arg3])
        return np.polyfit(x, y, 1)

    def fill_in_missing_data(self, old_df):
        """
        填补机器检测造成的遗漏数据错误填补成0的情况。
        """
        df = old_df.copy()
        for ion, series in df.iteritems():
            i = 0
            while i < len(series) - 2:
                correction_series = series[[df.index[i], df.index[i + 1], df.index[i + 2]]]  # 取数
                if correction_series.iloc[1] == 0 and correction_series.iloc[0] != 0 and correction_series.iloc[2] != 0:
                    series[df.index[i + 1]] = (correction_series.iloc[0] + correction_series.iloc[2]) / 2
                i += 1
        return df

    ############################
    def read_data(self, data):
        scan_index = np.array(data["scan_index"])
        intensity_values = np.array(data["intensity_values"])
        mass_values = np.array(data["mass_values"])
        scan_acquisition_time = np.array(data["scan_acquisition_time"])
        # 本段将原始数据转为行名rt，列名ion，值为响应的df，一个循环内的所有离子统一用循环开始时RT
        df = pd.DataFrame(index=scan_acquisition_time, columns=list(range(30, 401)))

        if len(scan_index) == len(intensity_values) == len(mass_values) == len(scan_acquisition_time):
            # 如果4 list长度相同，则用下列方法构建初始df
            r_mass_values = mass_values[::-1]
            w = 0
            while w <= len(mass_values) - 2:
                if r_mass_values[w + 1] - r_mass_values[w] < 0:
                    last_index = len(mass_values) - w - 1
                    break
                else:
                    w = w + 1
            # print("last_index = ", last_index)

            n = 0
            while n <= last_index:
                rt = scan_acquisition_time[n]
                p = n + 1
                if p > last_index:
                    break
                else:
                    while n < p <= last_index:
                        if mass_values[p + 1] - mass_values[p] > 0:
                            temp_ion = mass_values[n: p]
                            temp_intens = intensity_values[n: p]
                            for ion in temp_ion:
                                ion = int(ion)
                                m = np.where(temp_ion == ion)[0][0]
                                df.loc[rt, ion] = temp_intens[m]
                            n = p + 1
                            break
                        else:
                            p = p + 1

            # 下为将最后一个SIM循环添加到df中
            rt = scan_acquisition_time[last_index]
            temp_ion = mass_values[last_index:]
            temp_intens = intensity_values[last_index:]
            for ion in temp_ion:
                ion = int(ion)
                m = np.where(temp_ion == ion)[0][0]
                df.loc[rt, ion] = temp_intens[m]

        else:
            # 如果4 list长度不相同，则用下列方法构建初始df
            n = 0
            while n < len(scan_index) - 1:
                rt = scan_acquisition_time[n]
                index = scan_index[n]
                index_1 = scan_index[n + 1]
                temp_ion = mass_values[index:index_1]
                temp_intens = intensity_values[index:index_1]
                for ion in temp_ion:
                    ion = int(ion)
                    m = np.where(temp_ion == ion)[0][0]
                    df.loc[rt, ion] = temp_intens[m]
                n += 1

            rt = scan_acquisition_time[-1]
            index = scan_index[-1]
            temp_ion = mass_values[index:]
            temp_intens = intensity_values[index:]
            for ion in temp_ion:
                ion = int(ion)
                m = np.where(temp_ion == ion)[0][0]
                df.loc[rt, ion] = temp_intens[m]
        df = df.dropna(axis=1, how="all")
        # 去掉全为NaN的列
        df = df.dropna(axis=0, how="all")
        # 去掉全为NaN的行
        df = df.fillna(0)
        df = self.fill_in_missing_data(df)
        df.index = df.index.astype(float)
        return df

    def read_data_mzML(self, filename):
        run = pymzml.run.Reader(filename)
        # 创建一个空的DataFrame
        df = pd.DataFrame()

        # 遍历质谱谱图
        for spectrum in run:
            # 获取当前扫描的保留时间
            # print("*****", spectrum.peaks('centroided'))
            rt = (spectrum.scan_time_in_minutes()) * 60
            # 遍历每个峰，将峰的m/z值和强度保存到字典中
            mz_int_dict = {}
            for mz, intensity in spectrum.peaks('centroided'):
                mz_int_dict[int(mz)] = intensity

            # 将字典中的数据添加到DataFrame中
            row = pd.DataFrame(mz_int_dict, index=[rt])
            df = pd.concat([df, row], sort=True)

        df = df.dropna(axis=1, how="all")
        df = df.dropna(axis=0, how="all")
        df = df.fillna(0)
        df = self.fill_in_missing_data(df)
        df.index = df.index.astype(float)

        return df

    def export_msp(self, group_rt_s, peak_group_df):
        idx = list(peak_group_df.rt.values).index(group_rt_s)
        msp_df = peak_group_df.ions.values[idx]
        # 将DataFrame保存为以空格分隔的文本文件，并添加额外的文本信息
        with open(f'Group_{group_rt_s}.msp', 'w') as file:
            file.write(f"Name: Group_{group_rt_s}\n")
            file.write(f"Num peaks: {msp_df.shape[0]}\n")
            msp_df.to_csv(file, sep=' ', index=True, header=False, line_terminator='\n')

    def Main(self, total_result_pkl, peak_group_pkl, filename, smooth_value, peak_filt_value, group_peak_factor, msp,
             chooseinfo, RTRIinfo,
             page3none_match_weight, page3none_r_match_weight, page3none_group_weight, page3none_direct_weight,
             page3none_group_minimum_number_of_ions, page3none_sim_threshold, page3RT_search_wid,
             page3RT_match_weight, page3RT_r_match_weight, page3RT_group_weight, page3RT_direct_weight,
             page3RT_group_minimum_number_of_ions, page3RT_ri_participate, page3RT_sim_threshold, page3RT_window,
             page3RT_level_factor, page3RT_max_penalty, page3RT_no_info_penalty,
             page3RI_search_wid, page3RI_match_weight,
             page3RI_r_match_weight, page3RI_group_weight, page3RI_direct_weight, page3RI_group_minimum_number_of_ions,
             page3RI_RI_max, page3RI_ri_participate, page3RI_sim_threshold, page3RI_window,
             page3RI_window_scale, page3RI_level_factor, page3RI_max_penalty, page3RI_no_info_penalty,
             page3RI_inaccurate_ri_threshold, page3RI_inaccurate_ri_level_factor):
        # yuannote: RI最小值默认为0，不做输入窗口

        RI_min = 0
        peak_group_df = pd.DataFrame()  # liunote return需要
        qualitative_and_quantitative_analysis_result = pd.DataFrame()
        smooth_df = pd.DataFrame()
        # total_result.pkl是否有存在必要？该处需修改，后面再说
        # yuannote：此处是file-open project接口，只有total_result.pkl存在时，project才能读取
        # 有3个临时文件，total_result.pkl是总的，peak_group是到峰聚类为止的数据
        # qualitative_and_quantitative_analysis.pkl是定性后的数据

        total_result_file = total_result_pkl  # liunote
        total_result = None

        if os.path.exists(total_result_file):  # liunote
            with open(total_result_file, "rb") as f:
                total_result = pickle.load(f)

                peak_group_df = total_result[0][4]
                qualitative_and_quantitative_analysis_result = total_result[1]
                smooth_df = total_result[0][0]

        if total_result is None:
            ##################
            # yuannote：此处是Re-analysis——All processing的接口，如果点击All processing，那么就删除total_result.pkl、peak_group.pkl
            # 以及qualitative_and_quantitative_analysis.pkl，从这里重新分析，注意后续的参数也要重新读取
            #peak_group_filename = peak_group_pkl  # liunote
            peak_group_result = None

            # if os.path.exists(peak_group_filename):
            #     with open(peak_group_filename, "rb") as f:
            #         peak_group_result = pickle.load(f)
            #         smooth_df, smooth_df_final, decon_peak_dic, decon_data_df, peak_group_df, quant_result_df = peak_group_result

            if peak_group_pkl != '':
                smooth_df, smooth_df_final, decon_peak_dic, decon_data_df, peak_group_df, quant_result_df, df = peak_group_pkl
                peak_group_result = peak_group_pkl

            #if peak_group_result is None:
                # yuannote: 此处是file-opendata接口
                # 读取所输入文件的后缀名，或要求用户选择数据类型：如果是mzML，则用read_data_mzML函数
                # 如果是cdf，则用read_data函数
            else:
                # filename = 'tomato0418.mzML'
                if isinstance(filename, pd.DataFrame):
                    df = filename
                else:
                    if filename.endswith(".mzML"):  # liunote
                        df = self.read_data_mzML(filename)
                    elif filename.endswith(".cdf"):
                        df = self.read_data(filename)
                RT_max = max(df.index) / 60
                # yuannote: lwma里的5是向导页面1里的平滑系数
                smooth_df = self.lwma(df, smooth_value)  # liunote
                # yuannote: derivative里的10是指当ff等参数接近0时，用10替换，该参数默认10即可，不做输入接口
                smooth_df, noise_df = self.derivative(smooth_df, 10)
                ion_list = df.columns
                # ion_list = [41, 42, 43, 44, 45, 46, 47, 48]
                # del df
                # yuannote: find_peak里的10是是向导页面1里的峰过滤系数
                peak_dic, con_peak_dic, smooth_df_final = self.find_peak(ion_list, smooth_df, noise_df,
                                                                         peak_filt_value)  # todo:1 -> 10 #liunote
                decon_peak_dic, decon_data_df = self.decon(ion_list, smooth_df_final, con_peak_dic)
                peak_dic, decon_peak_dic = self.calculate_sv(smooth_df_final, peak_dic, decon_data_df, decon_peak_dic)
                # yuannote: group_peak里倒数第二个数字：0.5（bin_num)，对应的是向导页面1里的峰聚类灵敏度
                matched_wave, peak_group_df, quant_result_df = self.group_peak(peak_dic, decon_peak_dic, smooth_df,
                                                                               smooth_df_final,
                                                                               decon_data_df, 5,
                                                                               0.5, group_peak_factor, 0)  # liunote
                # yuannote：peak_group_df就是“软件任务栏下方的界面”中，“列1”框里要显示的数据
                # 其中：“只显示已定性group（打钩）”和“显示group的最少离子数"只改变在可视化界面显示的数据，不改变数据
                # 怎么改变在ppt里有，如果不清楚问我
                del peak_dic, noise_df

                peak_group_result = smooth_df, smooth_df_final, decon_peak_dic, decon_data_df, peak_group_df, quant_result_df, df
                # with open(peak_group_filename, "wb") as f:
                #     pickle.dump(peak_group_result, f)
            gc.collect()
            # yuannote：此处是Re-analysis——Identification的接口，如果点击Identification，
            # 那么就从读取已有的peak_group.pkl，删除已有的qualitative_and_quantitative_analysis.pkl和total_result.pkl
            # 重新进行定性分析，注意参数重新读取
            qualitative_and_quantitative_analysis_filename = "qualitative_and_quantitative_analysis.pkl"
            qualitative_and_quantitative_analysis_result = None

            if os.path.exists(qualitative_and_quantitative_analysis_filename):
                with open(qualitative_and_quantitative_analysis_filename, "rb") as f:
                    qualitative_and_quantitative_analysis_result = pickle.load(f)
            if qualitative_and_quantitative_analysis_result is None:
                # msp_file = r"./Remove_Duplicates.msp"  #liunote
                # yuannote: read_msp对应向导页面2的读取msp按钮
                meta, RI_df = self.read_msp(msp)
                # yuannote: RT_df_filename对应：向导页面2，基于保留信息分析数据选择RT时，点击读取保留信息按钮，读取RT信息
                RT_df_filename = ''  # liunote
                RI_data_filename = ''  # liunote
                if chooseinfo == "基于RT":  # liunote
                    RT_df_filename = RTRIinfo  # 用户输入的RT
                # yuannote: RI_data_filename：向导页面2，基于保留信息分析数据选择RI时，点击读取保留信息按钮，读取RI信息
                elif chooseinfo == "基于RI":  # liunote
                    RI_data_filename = RTRIinfo  # 用户输入的RI
                # species_pool_filename = r"02.csv"
                # species_pool = pd.DataFrame()
                if RT_df_filename != '':  # liunote
                    # yuannote：这个if表示向导页面2，基于保留信息分析数据选择RT时，读取了RT信息，进到了“RT-向导页面3”
                    # 如果向导页面2，基于保留信息分析数据选择RT时，未读取RT信息，则不能进行下一步，要求必须读取，ppt里有标注

                    RT_df = pd.read_csv(RT_df_filename, sep=",", index_col=0)
                    duplicated_index = RT_df[RT_df.index.duplicated()].index

                    error_df = RT_df.loc[duplicated_index, :].copy()

                    RT_df = RT_df[~RT_df.index.duplicated(keep='first')]
                    print("duplicated name in RT_list:")
                    print(error_df)

                    # 不排序的话find_nearest有问题
                    RT_df = RT_df.sort_values(by='RT', ascending=True)

                    # error_df.to_excel('./RT_error.xlsx')

                    # if os.path.exists(species_pool_filename):
                    #     species_pool = pd.read_csv(species_pool_filename, sep=",")

                    # yuannote:这里的参数都是“RT-向导页面3”中的参数，在后面做了注释

                    search_wid = page3RT_search_wid  # 搜库范围 #liunote
                    match_weight = page3RT_match_weight  # Match系数
                    r_match_weight = page3RT_r_match_weight  # R_Match系数
                    group_weight = page3RT_group_weight  # 峰组打分权重
                    direct_weight = page3RT_direct_weight  # 直接搜索打分权重
                    minimum_number_of_ions = page3RT_group_minimum_number_of_ions  # 参与定性group所含最少离子数
                    rt_participate = page3RT_ri_participate
                    # yuannote: RT-向导页面3中的RT参与打分（打钩），该选项打钩后，下面的参数才可设置
                    # 打钩即ri_participate = 1
                    sim_threshold = page3RT_sim_threshold
                    if rt_participate == True:  # liunote
                        window = page3RT_window  # 罚分窗口
                        level_factor = page3RT_level_factor  # 罚分系数
                        max_penalty = page3RT_max_penalty  # 最大罚分
                        no_info_penalty = page3RT_no_info_penalty  # 保留信息缺失罚分
                    if rt_participate == False:
                        # yuannote: RT-向导页面3中的RT参与打分（打钩），该选项不打钩时，即rt_participate == 0
                        # 不打钩时，罚分窗口、罚分系数、最大罚分、保留信息缺失罚分都是灰色的，不可更改
                        window = 0.3  # 罚分窗口
                        level_factor = 0.05  # 罚分系数
                        max_penalty = 0  # 最大罚分
                        no_info_penalty = 0  # 保留信息缺失罚分

                    qualitative_and_quantitative_analysis_result = self.peak_group_identification(peak_group_df,
                                                                                                  smooth_df,
                                                                                                  smooth_df_final, meta,
                                                                                                  RT_df,
                                                                                                  decon_peak_dic,
                                                                                                  decon_data_df, 1,
                                                                                                  match_weight,
                                                                                                  r_match_weight,
                                                                                                  search_wid,
                                                                                                  "RT", window, 5,
                                                                                                  level_factor,
                                                                                                  max_penalty,
                                                                                                  no_info_penalty,
                                                                                                  minimum_number_of_ions,
                                                                                                  sim_threshold,
                                                                                                  group_weight,
                                                                                                  direct_weight,
                                                                                                  0,
                                                                                                  0
                                                                                                  )
                    # yuannote: qualitative_and_quantitative_analysis_result里包含“软件任务栏下方的界面”页里“列2”和“列3”所要显示的信息

                elif RI_data_filename != '':
                    # yuannote：这个if表示向导页面2，基于保留信息分析数据选择RI时，读取了RI信息，进到了“RI-向导页面3”
                    # 如果向导页面2，基于保留信息分析数据选择RI时，未读取RI信息，则不能进行下一步，要求必须读取，ppt里有标注

                    standard_df = pd.read_csv(RI_data_filename, sep=",")
                    # yuannote: 这里对应RI-向导页面3-定性参数-最大保留指数
                    RI_max = page3RI_RI_max  # liunote

                    peak_group_df = self.peak_group_df_add_retention_infor(peak_group_df, "Kovats", standard_df, RI_min,
                                                                           RI_max)

                    duplicated_index = RI_df[RI_df.index.duplicated()].index

                    error_df = RI_df.loc[duplicated_index, :].copy()

                    RI_df = RI_df[~RI_df.index.duplicated(keep='first')]
                    print("duplicated name in RI_list:")
                    print(error_df)

                    # 不排序的话find_nearest有问题
                    RI_df = RI_df.sort_values(by='RI', ascending=True)

                    # error_df.to_excel('./RI_error.xlsx')

                    # if os.path.exists(species_pool_filename):
                    #     species_pool = pd.read_csv(species_pool_filename, sep=",")

                    # yuannote:这里的参数都是“RI-向导页面3”中的参数，在后面做了注释
                    search_wid = page3RI_search_wid  # 搜库范围 #liunote
                    match_weight = page3RI_match_weight  # Match系数
                    r_match_weight = page3RI_r_match_weight  # R Match系数
                    group_weight = page3RI_group_weight  # 峰组打分权重
                    direct_weight = page3RI_direct_weight  # 直接搜索打分权重
                    minimum_number_of_ions = page3RI_group_minimum_number_of_ions  # 参与定性group所含最少离子数
                    ri_participate = page3RI_ri_participate  # 01控制是否参与
                    # yuannote: RI-向导页面3中的RI参与打分（打钩），该选项打钩后，下面的参数才可设置
                    # 打钩即ri_participate = 1
                    sim_threshold = page3RI_sim_threshold
                    print("sim_threshod--------", sim_threshold)
                    print("type(sim_threshold---)", type(sim_threshold))
                    if ri_participate == True:  # liunote
                        window = page3RI_window  # 罚分窗口
                        window_scale = page3RI_window_scale  # RI窗口调整系数
                        level_factor = page3RI_level_factor  # 罚分系数
                        max_penalty = page3RI_max_penalty  # 最大罚分
                        no_info_penalty = page3RI_no_info_penalty  # 保留信息缺失罚分
                        inaccurate_ri_threshold = page3RI_inaccurate_ri_threshold
                        inaccurate_ri_level_factor = page3RI_inaccurate_ri_level_factor

                    if ri_participate == False:
                        # yuannote: RT-向导页面3中的RI参与打分（打钩），该选项不打钩时，即rt_participate == 0
                        # 不打钩时，罚分窗口、RI窗口调整系数、罚分系数、最大罚分、保留信息缺失罚分都是灰色的，不可更改
                        # 以下为固定参数，不可自定义
                        window = 100  # 罚分窗口
                        window_scale = 0  # RI窗口调整系数
                        level_factor = 0  # 罚分系数
                        max_penalty = 0  # 最大罚分
                        no_info_penalty = 0  # 保留信息缺失罚分
                        inaccurate_ri_threshold = 0
                        inaccurate_ri_level_factor = 0

                    qualitative_and_quantitative_analysis_result = self.peak_group_identification(peak_group_df,
                                                                                                  smooth_df,
                                                                                                  smooth_df_final, meta,
                                                                                                  RI_df,
                                                                                                  decon_peak_dic,
                                                                                                  decon_data_df, 1,
                                                                                                  match_weight,
                                                                                                  r_match_weight,
                                                                                                  search_wid,
                                                                                                  "RI", window,
                                                                                                  window_scale,
                                                                                                  level_factor,
                                                                                                  max_penalty,
                                                                                                  no_info_penalty,
                                                                                                  minimum_number_of_ions,
                                                                                                  sim_threshold,
                                                                                                  group_weight,
                                                                                                  direct_weight,
                                                                                                  inaccurate_ri_threshold,
                                                                                                  inaccurate_ri_level_factor
                                                                                                  )
                    # yuannote: qualitative_and_quantitative_analysis_result里包含“软件任务栏下方的界面”页里“列2”和“列3”所要显示的信息

                else:
                    # yuannote：这个if表示向导页面2，基于保留信息分析数据选择不使用时，进到了“不使用-向导页面3”

                    # yuannote:这里的参数都是“不使用-向导页面3”中的参数，在后面做了注释
                    match_weight = page3none_match_weight  # Match系数 #liunote
                    r_match_weight = page3none_r_match_weight  # R Match系数
                    group_weight = page3none_group_weight  # 峰组打分权重
                    direct_weight = page3none_direct_weight  # 直接搜索打分权重
                    minimum_number_of_ions = page3none_group_minimum_number_of_ions  # 参与定性group所含最少离子数
                    sim_threshold = page3none_sim_threshold
                    # if os.path.exists(species_pool_filename):
                    #     species_pool = pd.read_csv(species_pool_filename, sep=",")
                    qualitative_and_quantitative_analysis_result = self.peak_group_identification(peak_group_df,
                                                                                                  smooth_df,
                                                                                                  smooth_df_final, meta,
                                                                                                  RI_df,
                                                                                                  decon_peak_dic,
                                                                                                  decon_data_df, 1,
                                                                                                  match_weight,
                                                                                                  r_match_weight, 1,
                                                                                                  "None", 0, 0, 0, 0,
                                                                                                  0,
                                                                                                  minimum_number_of_ions,
                                                                                                  sim_threshold,
                                                                                                  group_weight,
                                                                                                  direct_weight,
                                                                                                  0,
                                                                                                  0
                                                                                                  )
                    # yuannote: Export results导出的csv就是qualitative_and_quantitative_analysis_result这个df
                    # yuannote: qualitative_and_quantitative_analysis_result里包含“软件任务栏下方的界面”页里“列2”和“列3”所要显示的信息

                # yuannote: peak_group_search指用group打分还是在group的rt上，将全部离子送去打分
                # 该功能调试用，就不设计输入窗口了，默认为1

                # with open(qualitative_and_quantitative_analysis_filename, "wb") as f:
                #     pickle.dump(qualitative_and_quantitative_analysis_result, f)

            total_result = (peak_group_result, qualitative_and_quantitative_analysis_result)

            # with open(total_result_file, "wb") as f:  # liunote
            #     pickle.dump(total_result, f)

        total_result_pkl = pickle.dumps(total_result)

        # smooth_df = smooth_df.iloc[:, 0:400]  # liunote: 只取前面的离子强度信息传到界面画图

        # peak_group_df = pd.read_csv(r"C:\Users\86724\Desktop\peak_group_df.csv", index_col=0)
        # qualitative_and_quantitative_analysis_result = pd.read_csv(r"C:\Users\86724\Desktop\qualitative_and_quantitative_analysis_result.csv", index_col=0)
        # smooth_df = pd.read_csv(r"C:\Users\86724\Desktop\smooth_df.csv", index_col=0)
        print("All_Done")
        # return peak_group_df, qualitative_and_quantitative_analysis_result, smooth_df
        return total_result_pkl


# a = DataAnalysis()
# total_result_pkl = 'total_result_pkl'
# peak_group_pkl = 'peak_group_pkl'
# filename = 'D:/work/GC方法/GC可视化/work/DataAnalysis/tomato0418.mzML'
# smooth_value = 5
# peak_filt_value = 10
# group_peak_factor = 1
# msp = 'D:/work/GC方法/GC可视化/work/DataAnalysis/Remove_Duplicates.msp'
# chooseinfo = "基于RT"
# RTRIinfo= 'D:/work/GC方法/GC可视化/work/DataAnalysis/New_RT_list.csv'
#
#
#
#
#
#
# page3none_match_weight = 0.3
# page3none_r_match_weight = 0.7
# page3none_group_weight = 0.2
# page3none_direct_weight = 0.8
# page3none_group_minimum_number_of_ions = 1
# page3none_sim_threshold = 0.3
#
# page3RT_search_wid = 1.5
# page3RT_match_weight = 0.3
# page3RT_r_match_weight = 0.7
# page3RT_group_weight = 0.2
# page3RT_direct_weight = 0.8
# page3RT_group_minimum_number_of_ions = 1
# page3RT_ri_participate = True
# page3RT_sim_threshold = 0.3
# page3RT_window = 0.30
# page3RT_level_factor = 0.05
# page3RT_max_penalty = 0.1
# page3RT_no_info_penalty = 0.05
#
#
#
# page3RI_search_wid = 150
# page3RI_match_weight = 0.7
# page3RI_r_match_weight = 0.3
# page3RI_group_weight = 0.2
# page3RI_direct_weight = 0.8
# page3RI_group_minimum_number_of_ions = 1
# page3RI_RI_max = 3000
# page3RI_ri_participate = True
# page3RI_sim_threshold = 0.4
# page3RI_window = 10
# page3RI_window_scale = 2
# page3RI_level_factor = 0.2
# page3RI_max_penalty = 0.4
# page3RI_no_info_penalty = 0.3
# page3RI_inaccurate_ri_threshold = 800
# page3RI_inaccurate_ri_level_factor = 0.01
#
# a.Main(total_result_pkl, peak_group_pkl, filename, smooth_value, peak_filt_value, group_peak_factor, msp, chooseinfo, RTRIinfo,
#              page3none_match_weight, page3none_r_match_weight, page3none_group_weight, page3none_direct_weight, page3none_group_minimum_number_of_ions, page3none_sim_threshold, page3RT_search_wid,
#        page3RT_match_weight, page3RT_r_match_weight, page3RT_group_weight, page3RT_direct_weight,
#        page3RT_group_minimum_number_of_ions, page3RT_ri_participate, page3RT_sim_threshold, page3RT_window,
#        page3RT_level_factor, page3RT_max_penalty, page3RT_no_info_penalty,
#        page3RI_search_wid, page3RI_match_weight,
#        page3RI_r_match_weight, page3RI_group_weight, page3RI_direct_weight, page3RI_group_minimum_number_of_ions,
#        page3RI_RI_max, page3RI_ri_participate, page3RI_sim_threshold, page3RI_window,
#        page3RI_window_scale, page3RI_level_factor, page3RI_max_penalty, page3RI_no_info_penalty,
#        page3RI_inaccurate_ri_threshold, page3RI_inaccurate_ri_level_factor
#        )

