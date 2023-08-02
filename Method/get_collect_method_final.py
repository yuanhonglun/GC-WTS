import pandas as pd
import numpy as np
import re
import os


class GetMethod():

    def replace_name(self, list_1):
        replacements = {".alpha.": "alpha", ".beta.": "beta", ".gamma.": "gamma", ".delta.": "delta",
                        ".omega.": "omega", ".tau.": "tau"}
        for i in range(len(list_1)):
            for old_str, new_str in replacements.items():
                list_1[i] = list_1[i].replace(old_str, new_str)

        return list_1

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

    def read_msp(self, msp_file):
        msp_file = open(msp_file, "r")
        list_1 = msp_file.readlines()
        # print(list_1)
        new_list = [item.replace('NAME:', 'Name:') for item in list_1]
        list_1 = [item.replace("Num peaks:", "Num Peaks:") for item in new_list]
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

        return meta

    def dot_product_distance(self, p, q):

        if (np.sum(p)) == 0 or (np.sum(q)) == 0:
            score = 0
        else:
            score = np.power(np.sum(q * p), 2) / \
                    (np.sum(np.power(q, 2)) * np.sum(np.power(p, 2)))
        return score

    def weighted_dot_product_distance(self, compare_df, fr_factor):

        # print("debug compare_df = ")
        # print(compare_df)
        m_q = pd.Series(compare_df.index)
        m_q = m_q.astype(float)
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
        if m >= fr_factor:
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

    def calculate_similarity(self, target_name, df, n, fr_factor):
        result_df = pd.DataFrame(columns=["Score"])
        first_col = df.loc[target_name]
        for compound in df.index.values:
            if compound != target_name:
                second_col = df.loc[compound]
                compare_df = pd.concat([first_col, second_col], axis=1)
                compare_df = compare_df.astype(float)
                # print(compound)
                # print(compare_df)
                score = self.weighted_dot_product_distance(compare_df, fr_factor)
                result_df.loc[compound, "Score"] = score

        return result_df

    def calculate_average_score_and_difference_count(self, targeted_compound,
                                                     ion_combination,
                                                     df, similarity_threshold,
                                                     fr_factor,
                                                     n):
        # 此处输入的df是已将丰度40以下的响应归零，且去掉目标物响应为0列的df
        # print("debug ion_combination = ", ion_combination)
        difference_count_df = pd.DataFrame(columns=["Diff_Count", "Similar_Compound_Ave_Score"])
        for ions in ion_combination:
            # print(ions)
            # ions = [str(i) for i in ions]
            temp_df_1 = df[ions]
            # temp_df_1 = temp_df_1.astype(int)
            # 根据x组合，切片dataframe
            result_df_2 = self.calculate_similarity(targeted_compound, temp_df_1, n, fr_factor)
            # yuannote: fr_factor，选离子循环时，5个或以上的共有离子才计算ss得分的fr部分，防止离子过少时fr大幅降低ss得分

            result_df_3 = result_df_2[(result_df_2["Score"] < similarity_threshold)]
            # result_df_3是将result_df_2中相似性小于0.85的提出来，作为一个新df
            count = len(result_df_3)
            # todo: 此处可能需修改
            difference_count_df.loc[str(ions), "Diff_Count"] = count

            result_df_4 = result_df_2[(result_df_2["Score"] >= similarity_threshold)]
            # print("尚未分清物质")
            # print(result_df_4)
            # result_df_4是将result_df_2中相似性大于等于0.85的提出来，作为一个新df
            # result_df_4的行名是现阶段尚未区分的物质名
            if result_df_4.shape[0] > 0:
                ave_score_1 = np.average(result_df_4, axis=0)[0]
            else:
                ave_score_1 = 1
            # todo: 此处可能需修改
            difference_count_df.loc[str(ions), "Similar_Compound_Ave_Score"] = ave_score_1
            # jf：极端值是否影响平均值

        difference_count_df.sort_values(by="Diff_Count", inplace=True, ascending=False)
        # difference_count_df.to_excel('C:/HNU/课题/高通量GC-MS方法/2.0/数据库+质谱数据/特异离子提取/results/{}.xlsx'.format("difference_count_" + ions), index=True)
        # ave_score_df.to_excel('C:/HNU/课题/高通量GC-MS方法/2.0/数据库+质谱数据/特异离子提取/results/{}.xlsx'.format("ave_score_" + ions), index=True)
        # 排序

        return difference_count_df

    def calculate_combination_score(self, combination_df, targeted_compound, temp_df, prefer_mz_threshold):
        for index, row in combination_df.iterrows():
            ion_list = re.findall("\d+\.?\d*", index)
            ion_list = list(map(int, ion_list))
            # print("debug_ion_combination = ", ion_list)
            new_temp_df = temp_df.loc[str(targeted_compound), ion_list].to_frame()
            new_temp_df["ion"] = new_temp_df.index.tolist()
            new_temp_df["ion"] = new_temp_df["ion"].astype('int')
            new_temp_df["ion"] = np.where(new_temp_df["ion"] < prefer_mz_threshold, 1, new_temp_df["ion"])
            new_temp_df["score"] = (pow(new_temp_df["ion"], 3)) * (pow(new_temp_df[str(targeted_compound)], 0.5))
            combination_df.loc[index, "com_score"] = new_temp_df["score"].sum()

        return combination_df

    def calculate_solo_compound_combination_score(self, matrix_1, prefer_mz_threshold):
        matrix_1['ion'] = matrix_1['ion'].apply(lambda x: 1 if x < prefer_mz_threshold else x)
        # matrix_1.to_excel(r".\matrix_1_阈值校正后.xlsx", index=True)
        matrix_1['com_score'] = matrix_1.apply(lambda row: pow(row.iloc[0], 0.5) * pow(row.iloc[1], 3), axis=1)
        matrix_1 = matrix_1.sort_values(by='com_score', ascending=False)
        # matrix_1.to_excel(r".\matrix_1_计算打分后.xlsx", index=True)

        return matrix_1

    def validateTitle(self, title):
        rstr = r"[\/\\\:\*\?\"\<\>\|]"  # '/ \ : * ? " < > |'
        new_title = re.sub(rstr, "_", title)  # 替换为下划线
        return new_title

    def group_rows(self, row):
        return '_'.join(row.astype(str))

    def replace_1(self, x):
        if 1 in x.values:
            x[:] = 1
        return x

    # YuanNote: 这里输入3个文件，msp库、RT列表，如果set_name为1，则需输入物质名列表，否则不输入

    def Main(self, msp_path, rt_data_path, set_name_list, name_list_path, mz_min, mz_max, outpath, rt_window,
             min_ion_intensity_percent, min_ion_num, prefer_mz_threshold, similarity_threshold, fr_factor,
             retention_time_max, solvent_delay, sim_sig_max, min_dwell_time,
             min_point_per_s, min_point_per_s_limit, convert_to_ag_method):

        msp = msp_path
        RT_data = pd.read_excel(rt_data_path, header=0, index_col=0)
        error_df = pd.DataFrame(columns=["error"])
        meta_1 = self.read_msp(msp)

        matrix_file = "matrix.xlsx"
        matrix = None
        if os.path.exists(matrix_file):
            matrix = pd.read_excel(r"./matrix.xlsx", header=0, index_col=0)

        if matrix is None:
            matrix = pd.DataFrame()
            for name, dic in meta_1.items():
                b = {name: dic}
                y = pd.DataFrame.from_dict(b).T
                matrix = pd.concat([matrix, y], axis=0, join='outer').fillna(0)

            for col_name in list(matrix):
                if type(col_name) not in [float, int, np.float64, np.int64] and col_name.isdigit() == False:
                    print(col_name)
                    matrix.drop(columns=col_name, inplace=True)
                    # 删除不是数字的列

            for ion in list(matrix):
                if int(ion) < mz_min or int(ion) > mz_max:
                    matrix.drop(columns=ion, inplace=(True))

            matrix.to_excel(r".\matrix.xlsx", index=True)

        matrix_name = matrix.index.tolist()

        duplicated_index = RT_data.index[RT_data.index.duplicated()]
        for index in duplicated_index:
            error_df.loc[index, "error"] = "duplicated_RT"
        RT_data.drop(duplicated_index, axis=0, inplace=True)
        # 找RT中重复，并导入到error_df中
        for index_1, row in RT_data.iterrows():

            if index_1 not in matrix_name:
                print(index_1, "in RT_list, but not in MSP library")
                error_df.loc[index_1, "error"] = "In RT list, but not in MSP library"
                RT_data.drop(index=index_1, inplace=True)
                # 从RT列表里去掉MSP库里没有的物质，并append到error_list上

            elif type(row[0]) not in [float, int, np.float64, np.int64] or row[0] == np.nan:
                print(type(row[0]))
                print("Error in RT format of ", index_1)
                error_df.loc[index_1, "error"] = "RT format error"
                RT_data.drop(index=index_1, inplace=True)

        RT_data = RT_data.sort_values(by='RT')

        # yuannote: 这里选择是否输入name_list
        if set_name_list:
            f1 = open(name_list_path)
            f2 = list(f1)
            compound_list = []
            for name in f2:
                name = name.strip()
                # 去掉物质名后面的/n
                name = name.strip('"')
                # 去掉物质两边的引号
                compound_list.append(name)
            # print(compound_list)
            compound_list = list(set(compound_list))
            for i in compound_list:
                if i not in RT_data.index.values:
                    print(i, "in name list, but not in RT list")
                    error_df.loc[i, "error"] = "in name list, but not in RT list"
                    compound_list.remove(i)

        else:
            compound_list = RT_data.index.values.tolist()

        # yuannote: 这里输出第一个结果
        error_df.to_excel(outpath + "/{}.xlsx".format("input_data_error_info"), index=True)

        # yuannote: 这里要求用户输入保留时间窗口，默认0.5
        # rt_window = 0.5
        nearby_compound_dic = {}
        for name in compound_list:
            if name in RT_data.index.values.tolist():
                rt = RT_data.at[name, "RT"]
                nearby_compound_dic[name] = RT_data[(RT_data.iloc[:, 0] >= rt - rt_window) &
                                                    (RT_data.iloc[:, 0] <= rt + rt_window)].index.tolist()
        # print("nearby_compound_dic = ", nearby_compound_dic)

        combination_result_df_file = "combination_results.xlsx"
        combination_result_df = None

        if os.path.exists(combination_result_df_file):
            combination_result_df = pd.read_excel(r"./combination_results.xlsx", header=0, index_col=0)

        if combination_result_df is None:
            combination_result_df = pd.DataFrame(
                columns=["RT", "Ion_Combination", "Note", "Similar_Compound_List", "SCL_Note"])

            # yuannote: 这里让用户输入选离子时的最小响应和最少离子数，注意最小响应输入的是7%，那么进来后
            # 变成70，只要求输入min_ion_intensity_percent和min_ion_num
            # min_ion_intensity_percent = 5
            min_ion_intensity = min_ion_intensity_percent * 10
            # min_ion_num = 2
            # prefer_mz_threshold = 60
            # similarity_threshold = 0.85
            # fr_factor = 2

            for targeted_compound, nearby_compound_list in nearby_compound_dic.items():
                combination_result_df.loc[targeted_compound, "RT"] = RT_data.loc[targeted_compound, "RT"]
                print("*" * 30)
                print(targeted_compound)
                title = self.validateTitle(str(targeted_compound))
                # 替换物质名中不能作为文件名的字符
                solo_list = []
                error_list_2 = []
                if nearby_compound_list == [targeted_compound]:
                    scl = []
                    combination_result_df.loc[targeted_compound, "Similar_Compound_List"] = scl
                    combination_result_df.loc[targeted_compound, "SCL_Note"] = "无临近物质"
                    print("similar_compound_list")
                    print(scl)

                    matrix_1 = matrix.loc[targeted_compound].to_frame()
                    matrix_1["ion"] = matrix_1.index.tolist()
                    matrix_1["ion"] = matrix_1["ion"].astype(int)
                    matrix_1[targeted_compound] = matrix_1[targeted_compound].astype(float)
                    matrix_1[targeted_compound] = np.where(matrix_1[targeted_compound] < min_ion_intensity, 0,
                                                           matrix_1[targeted_compound])
                    matrix_1 = matrix_1.loc[matrix_1[targeted_compound] > 0, :]

                    if matrix_1.shape[0] < 2:
                        print("final_combination = NA")
                        combination_result_df.loc[targeted_compound, "Ion_Combination"] = "NA"
                        combination_result_df.loc[targeted_compound, "Note"] = "无临近物质；经最小离子响应筛选后，可用离子数量小于2，不参与筛选"
                    else:
                        # matrix_1.to_excel(r".\matrix_1_原始.xlsx", index=True)
                        matrix_1 = self.calculate_solo_compound_combination_score(matrix_1, prefer_mz_threshold)

                        if matrix_1.shape[0] <= min_ion_num:
                            combination_list = matrix_1.index.values.tolist()
                        else:
                            combination_list = matrix_1.iloc[0: min_ion_num, :].index.values.tolist()

                        # 下为原solo物质离子选择算法，根据solo_compound_prefer_mz_threshold，选大于mz阈值响应最高的离子，不足用其他离子补足
                        # matrix_1_1 = matrix_1.loc[matrix_1.loc[:,"ion"] >= solo_compound_prefer_mz_threshold,:]
                        # matrix_1_2 = matrix_1.loc[matrix_1.loc[:,"ion"] < solo_compound_prefer_mz_threshold,:]
                        # matrix_1_1 = matrix_1_1.sort_values(by = str(targeted_compound), ascending=False, inplace=False)
                        # matrix_1_2 = matrix_1_2.sort_values(by = str(targeted_compound), ascending=False, inplace=False)
                        # matrix_1_3 = pd.concat((matrix_1_1, matrix_1_2))
                        # if matrix_1_3.shape[0] < min_ion_num :
                        #     combination_list = matrix_1_3.iloc[:, 1].tolist()
                        # else:
                        #     combination_list = matrix_1_3.iloc[0:min_ion_num, 1].tolist()
                        # print("debug_solo_compound_combination_list = ", combination_list)
                        combination_result_df.loc[str(targeted_compound), "Ion_Combination"] = combination_list
                        print("final_combination = ", combination_list)
                        solo_list.append(targeted_compound)
                        # 如果临近物质列表为空，那么就取该物质mz80以上的离子中，响应最强的3个，若没有或不够，用80以下最强响应离子补足
                else:
                    # print("nearby_compound_list = ", nearby_compound_list)
                    temp_df = matrix.loc[nearby_compound_list]
                    # 根据行名列表提取出新dataframe
                    temp_df = temp_df.astype(float)
                    # temp_df.to_excel(r".\原始temp_df.xlsx", index=True)
                    temp_df.loc[targeted_compound, :] = np.where(temp_df.loc[targeted_compound, :] < min_ion_intensity,
                                                                 0,
                                                                 temp_df.loc[targeted_compound, :])
                    # temp_df.to_excel(r".\去最小值temp_df.xlsx", index=True)
                    temp_df = temp_df.loc[:, temp_df.loc[targeted_compound, :] > 0]
                    # temp_df.to_excel(r".\去零temp_df.xlsx", index=True)
                    # 目标物中响应小于min_ion_intensity的离子，响应改为0，后去掉目标物中响应为0的离子，目的是选离子的时候不要低响应离子

                    if temp_df.shape[1] < 2:
                        print("final_combination = NA")
                        combination_result_df.loc[targeted_compound, "Ion_Combination"] = "NA"
                        combination_result_df.loc[targeted_compound, "Note"] = "经最小离子响应筛选后，可用离子数量小于2，不参与筛选"
                        # 这里是指经过最低响应过滤后，目标物剩余离子已小于2个，这种物质不参与选离子
                    else:
                        # temp_df = temp_df.replace(0, 1)
                        # temp_df.to_excel('C:/HNU/课题/高通量GC-MS方法/2.0/数据库+质谱数据/特异离子提取/results/temp_df.xlsx')
                        # 上面注释掉的是将其他物质中读数为0的离子替换成1，是否这样做跟后续的peak group
                        # 功能有关，详见思路0814.docx
                        similar_compound_list = []
                        # print("temp_df = ", temp_df)
                        result_df_1 = self.calculate_similarity(targeted_compound, temp_df, -1, fr_factor)
                        # yuannote: fr_factor，计算初始相似性时，2个或以上共有离子就计算fr得分，与非靶一致

                        # result_df_1.to_excel(r".\初始相似性判断.xlsx", index=True)
                        # print("初始相似性判断：")
                        # print(result_df_1)
                        # result_df_1.to_excel('C:/HNU/课题/高通量GC-MS方法/2.0/数据库+质谱数据/特异离子提取/results/{}.xlsx'.format(title + "_初始相似度"), index=True)
                        for index, row in result_df_1.iterrows():
                            if float(row) >= similarity_threshold:
                                similar_compound_list.append(index)
                                # 离子完整时，找出不能区分的物质(相似性>=0.85)，输入到similar_compound_list中备用

                        # print("*" * 30)
                        # print(targeted_compound)
                        # print("完整离子相似性计算")
                        # print(result_df_1)
                        # print("similar_compound_list")
                        # print(similar_compound_list)
                        combination_result_df.loc[targeted_compound, "Similar_Compound_List"] = similar_compound_list
                        temp_df.drop(index=similar_compound_list, inplace=True)
                        # 不能区分的物质，要从df里删除该行

                        if temp_df.shape[0] == 1:
                            temp_name = ((temp_df.index.values).tolist())[0]
                            if temp_name == targeted_compound:
                                matrix_1 = matrix.loc[targeted_compound].to_frame()
                                matrix_1["ion"] = matrix_1.index.tolist()
                                matrix_1["ion"] = matrix_1["ion"].astype(int)
                                matrix_1[targeted_compound] = matrix_1[targeted_compound].astype(float)
                                matrix_1[targeted_compound] = np.where(matrix_1[targeted_compound] < min_ion_intensity,
                                                                       0,
                                                                       matrix_1[targeted_compound])
                                matrix_1 = matrix_1.loc[matrix_1[targeted_compound] > 0, :]

                                if matrix_1.shape[0] < 2:
                                    print("final_combination = NA")
                                    combination_result_df.loc[targeted_compound, "Ion_Combination"] = "NA"
                                    combination_result_df.loc[targeted_compound, "Note"] = "经最小离子响应筛选后，可用离子数量小于2，不参与筛选"
                                else:
                                    # matrix_1.to_excel(r".\matrix_1_原始.xlsx", index=True)
                                    matrix_1 = self.calculate_solo_compound_combination_score(matrix_1,
                                                                                              prefer_mz_threshold)

                                    if matrix_1.shape[0] <= min_ion_num:
                                        combination_list = matrix_1.index.values.tolist()
                                    else:
                                        combination_list = matrix_1.iloc[0: min_ion_num, :].index.values.tolist()

                                    # matrix_1_1 = matrix_1.loc[matrix_1.loc[:,"ion"] >= solo_compound_prefer_mz_threshold,:]
                                    # matrix_1_2 = matrix_1.loc[matrix_1.loc[:,"ion"] < solo_compound_prefer_mz_threshold,:]
                                    # matrix_1_1 = matrix_1_1.sort_values(by = str(targeted_compound), ascending=False, inplace=False)
                                    # matrix_1_2 = matrix_1_2.sort_values(by = str(targeted_compound), ascending=False, inplace=False)
                                    # matrix_1_3 = pd.concat((matrix_1_1, matrix_1_2))
                                    # if matrix_1_3.shape[0] < min_ion_num or min_ion_num > matrix_1_3.shape[0]:
                                    #     combination_list = matrix_1_3.iloc[:, 1].tolist()
                                    # else:
                                    #     combination_list = matrix_1_3.iloc[0:min_ion_num, 1].tolist()
                                    # print("debug_solo_compound_combination_list = ", combination_list)
                                    combination_result_df.loc[
                                        str(targeted_compound), "Ion_Combination"] = combination_list
                                    print("final_combination = ", combination_list)
                                    solo_list.append(targeted_compound)
                                # 如果去掉similar_compound_list只剩自己，那么就取该物质mz80以上的离子中，响应最强的两个，若没有或不够，用80以下最强响应离子补足
                        else:

                            # 开始随机取2个离子变成组合
                            col_name = list(temp_df.columns)
                            col_name = list(map(int, col_name))
                            # print(col_name)

                            ########初始组合2离子算法########
                            # combination = list(itertools.product(col_name, repeat=2))
                            # # repeat = 2，即在col_name这一列表中两两组合，生成combination这一list，里面是('15', '18')这种元组(tuple)
                            # # print(type(combination))
                            # for x in combination:
                            #     if len(set(x)) != len(x):
                            #         combination.remove(x)
                            #         # 去掉了combination中元素重复的元组，例如('15', '15')
                            # new_com = []
                            # for x in combination:
                            #     x = list(x)
                            #     x.sort()
                            #     if x not in new_com:
                            #         new_com.append(x)
                            # # 去掉了(15, 17),(17, 15)这种重复，只保留一个
                            ########初始组合2离子算法########

                            ########初始组合1离子算法########
                            new_com = [[x] for x in col_name]
                            # print("new_com = ", new_com)
                            ########初始组合1离子算法########

                            difference_count_df_1 = self.calculate_average_score_and_difference_count(targeted_compound,
                                                                                                      new_com,
                                                                                                      temp_df,
                                                                                                      similarity_threshold,
                                                                                                      fr_factor,
                                                                                                      0)

                            # 开始取初始组合
                            combination_df = difference_count_df_1[difference_count_df_1["Diff_Count"] >=
                                                                   difference_count_df_1.iat[0, 0]]
                            if combination_df.shape[0] > 5:
                                combination_df = combination_df.sort_values(by="Similar_Compound_Ave_Score",
                                                                            inplace=False, ascending=True)
                                combination_df = combination_df[combination_df["Similar_Compound_Ave_Score"] >=
                                                                combination_df.iat[0, 1]]
                                # 取scas最低分的组合，可能多个
                                if combination_df.shape[0] > 5:
                                    combination_df = self.calculate_combination_score(combination_df, targeted_compound,
                                                                                      temp_df, prefer_mz_threshold)
                                    combination_df = combination_df.sort_values(by="com_score", inplace=False,
                                                                                ascending=False)
                                    combination_df = combination_df[:5]
                                else:
                                    combination_df = combination_df
                            else:
                                combination_df = combination_df
                            # 在ave_score最低的组合中，根据mass方×intensity开根打分，取前五，进入第一轮

                            # print("initial_com_df = ", combination_df)
                            # combination_df.to_excel('C:/HNU/课题/高通量GC-MS方法/2.0/数据库+质谱数据/特异离子提取/results/{}.xlsx'.format("top_combination_" + "round_0"), index=True)
                            ion_list = list(temp_df)
                            # ion_list = list(map(int, ion_list))
                            # 读取temp_df列名并转换为int类型
                            combination_array = combination_df.index.values
                            n = 0
                            ion_num = 1
                            # yuannote: 这里的ion_num跟初始取1离子还是2离子有关，改成一致的

                            flag = True

                            while True:
                                if int((combination_df.max())[0]) >= int(
                                        temp_df.shape[0] - 1) and ion_num >= min_ion_num:
                                    break
                                elif flag == False:
                                    break
                                else:
                                    # print(int((combination_df.max())[0]))
                                    n = n + 1
                                    total_list = []
                                    new_total = []
                                    # print("combination_array = ", combination_array)
                                    # print(type(combination_array))
                                    for ion_combination in combination_array:
                                        # print("ion_combination = ", ion_combination)
                                        # print(type(ion_combination))
                                        ion_combination_list = re.findall("\d+\.?\d*", ion_combination)
                                        ion_combination_list = list(map(int, ion_combination_list))
                                        # 将ion_combination(str类型)中的数字提取出来，转为list，list里元素转为int
                                        # print("ion_combination = ", ion_combination)
                                        # print("ion_list = ", ion_list)
                                        # print("ion_combination_list = ", ion_combination_list)
                                        candidate_list = [i for i in ion_list if i not in ion_combination_list]
                                        # 候选离子list = 总离子list-离子组合list
                                        # print("candidate_list = ", candidate_list)
                                        if candidate_list == []:
                                            if int((combination_df.max())[0]) >= int(temp_df.shape[0] - 1):
                                                print("离子数小于最小离子数，已实现分离")
                                                if combination_df.shape[0] > 1:
                                                    combination_df = combination_df.sort_values(
                                                        by="Similar_Compound_Ave_Score", inplace=False, ascending=True)
                                                    # 此处会报“A value is trying to be set on a copy of a slice from a DataFrame”的错误，将inplace从True
                                                    # 改为False，停止报错，即将变量更改后，赋予到原变量中
                                                    combination_df = combination_df[
                                                        combination_df["Similar_Compound_Ave_Score"] >=
                                                        combination_df.iat[0, 1]]
                                                    if combination_df.shape[0] > 1:
                                                        combination_df = self.calculate_combination_score(
                                                            combination_df,
                                                            targeted_compound,
                                                            temp_df,
                                                            prefer_mz_threshold)
                                                        combination_df = combination_df.sort_values(by="com_score",
                                                                                                    inplace=False,
                                                                                                    ascending=False)
                                                        combination_df = combination_df[:1]
                                                    else:
                                                        combination_df = combination_df
                                                else:
                                                    combination_df = combination_df
                                                combination_array = combination_df.index.values
                                                combination_result_df.loc[str(targeted_compound), "Ion_Combination"] = \
                                                    combination_array[0]
                                                combination_result_df.loc[
                                                    str(targeted_compound), "Note"] = "离子数小于最小离子数，已实现分离"
                                                print("final_combination = ", combination_array[0])
                                                flag = False

                                                break

                                            else:
                                                print("无离子可用，未实现分离")
                                                combination_result_df.loc[
                                                    str(targeted_compound), "Ion_Combination"] = "NA"
                                                combination_result_df.loc[
                                                    str(targeted_compound), "Note"] = "无离子可用，未实现分离"
                                                print("final_combination = NA")
                                                flag = False
                                                # 已用了所有离子，仍未分清物质的情况

                                                break

                                        elif flag == True:

                                            for candidate in candidate_list:
                                                # print("Before = ", ion_combination_list)
                                                temp_ion_combination_list = ion_combination_list.copy()
                                                temp_ion_combination_list.append(candidate)
                                                # print("After = ", ion_combination_list)
                                                total_list.append(temp_ion_combination_list)
                                                temp_ion_combination_list = []
                                            for i in total_list:
                                                i = list(map(int, i))
                                                # print(i)
                                                i.sort()
                                                if i not in new_total:
                                                    new_total.append(i)
                                                # 生成下一轮离子组合后，进行去重
                                        # print("new_total = ", new_total)

                                    if flag == True:

                                        difference_count_df_2 = self.calculate_average_score_and_difference_count(
                                            targeted_compound,
                                            new_total,
                                            temp_df,
                                            similarity_threshold,
                                            fr_factor,
                                            n)

                                        if len(difference_count_df_2) > 0:

                                            combination_df = difference_count_df_2[
                                                difference_count_df_2["Diff_Count"] >=
                                                difference_count_df_2.iat[0, 0]]
                                        else:
                                            print("bug: difference_count_df为空")
                                            combination_result_df.loc[str(targeted_compound), "Ion_Combination"] = "NA"
                                            combination_result_df.loc[
                                                str(targeted_compound), "Note"] = "bug: difference_count_df为空"
                                            print("final_combination = NA")
                                            flag = False

                                            break

                                        if combination_df.shape[0] > 1:
                                            combination_df = combination_df.sort_values(by="Similar_Compound_Ave_Score",
                                                                                        inplace=False, ascending=True)
                                            # 此处会报“A value is trying to be set on a copy of a slice from a DataFrame”的错误，将inplace从True
                                            # 改为False，停止报错，即将变量更改后，赋予到原变量中
                                            combination_df = combination_df[
                                                combination_df["Similar_Compound_Ave_Score"] >=
                                                combination_df.iat[0, 1]]
                                            if combination_df.shape[0] > 1:
                                                combination_df = self.calculate_combination_score(combination_df,
                                                                                                  targeted_compound,
                                                                                                  temp_df,
                                                                                                  prefer_mz_threshold)
                                                combination_df = combination_df.sort_values(by="com_score",
                                                                                            inplace=False,
                                                                                            ascending=False)
                                                combination_df = combination_df[:1]
                                            else:
                                                combination_df = combination_df
                                        else:
                                            combination_df = combination_df
                                            # 每轮的combination计算完，先判断是否能完全分离且离子数量大于最小离子数，是则停，否则计算每个组合与剩余物质相似性，取相似性最低组合，
                                        # combination_df.to_excel('C:/HNU/课题/高通量GC-MS方法/2.0/数据库+质谱数据/特异离子提取/results/{}.xlsx'.format("top_combination_" + "round_" + str(n)), index=True)
                                        combination_array = combination_df.index.values
                                        ion_num = ion_num + 1
                                        # print(ion_num)

                                    else:

                                        break
                                # print("第", n, "轮")
                                # print(combination_df)
                            if flag == True:
                                combination_result_df.loc[str(targeted_compound), "Ion_Combination"] = \
                                    combination_array[0]
                                print("final_combination = ", combination_array[0])

            # yuannote: 这里输出第二个结果
            combination_result_df.to_excel(outpath + "/combination_results.xlsx", index=True)

        error_df = pd.DataFrame(columns=["Name", "Error"])
        name_list_total = []
        num = []
        name_list = combination_result_df.index.values.tolist()
        RT_list_total = []
        for name in name_list:
            if name in RT_data.index.values.tolist():
                if type(combination_result_df.loc[name, "Ion_Combination"]) == str:
                    ion_str = combination_result_df.loc[name, "Ion_Combination"]
                    ion_list = re.findall("\d+\.?\d*", ion_str)
                    ion_list = list(map(int, ion_list))
                    for x in range(0, len(ion_list)):
                        name_list_total.append(name)
                        RT_list_total.append(RT_data.loc[name, "RT"])
                        num.append(ion_list[x])
                elif type(combination_result_df.loc[name, "Ion_Combination"]) == list:
                    ion_list = [int(x) for x in combination_result_df.loc[name, "Ion_Combination"]]
                    for x in range(0, len(ion_list)):
                        name_list_total.append(name)
                        RT_list_total.append(RT_data.loc[name, "RT"])
                        num.append(ion_list[x])
                else:
                    print(name, "Error in Ion_Combination format")
                    error_df.loc[len(error_df)] = [name, "离子组合格式有误"]
            else:
                error_df.loc[len(error_df)] = [name, "物质不在RT列表中"]

        length = len(name_list_total)

        data = {'RT': RT_list_total, 'ion': num}

        ion_rt = pd.DataFrame(data, index=name_list_total)

        # yuannote: 这里要求用户输入保留时间最大值
        # retention_time_max = 68.8

        ion_rt.sort_values(by="RT", inplace=True, ascending=True)
        for idx, row in ion_rt.iterrows():
            if row["RT"] > retention_time_max:
                ion_rt.loc[idx, "RT"] = retention_time_max

        # yuannote: 这里输出结果
        ion_rt.to_excel(outpath + "/ion_rt_data.xlsx", index=True)

        # yuannote: 这里输出第三个结果
        error_df.to_excel(outpath + "/converter_error_info.xlsx", index=False)

        ##############进行SIM分段###############

        rt_index = [i * 0.5 / 60 for i in range(0, int(retention_time_max * 120) + 1, 1)]

        # solvent_delay = 0
        # yuannote: 这里输入溶剂延迟，RT上小于溶剂延迟的离子不检测

        df = pd.DataFrame(index=rt_index, columns=[i for i in range(mz_min, mz_max + 1)])

        for i, row in ion_rt.iterrows():
            df.loc[(df.index > row[0] - rt_window) & (df.index < row[0] + rt_window), row[1]] = 1

        df = df[df.index > solvent_delay]

        df["sum"] = df.sum(axis=1)

        # 将各行的值转换为字符串并分组
        df['group'] = df.apply(self.group_rows, axis=1)

        df["group_id"] = ""
        df.iloc[0, -1] = 1
        n = 1
        m = 2
        while n <= len(df) - 1:
            if df.iloc[n, -2] == df.iloc[n - 1, -2]:
                df.iloc[n, -1] = df.iloc[n - 1, -1]
                n = n + 1
            else:
                df.iloc[n, -1] = m
                n = n + 1
                m = m + 1

        group_list = df['group_id'].unique().tolist()

        for i in group_list:
            if i in df['group_id'].values.tolist():
                first_row_sum = df.loc[df['group_id'] == i].iloc[0, -3]
                if first_row_sum == 0:
                    if i - 1 and i + 1 in group_list:
                        if df.loc[df['group_id'] == i - 1].iloc[0, -3] > df.loc[
                            df['group_id'] == i + 1].iloc[0, -3]:
                            mask = (df['group_id'] == i) | (df['group_id'] == i + 1)
                            selected_rows = df.loc[mask]
                            selected_rows = selected_rows.iloc[:, :-3].apply(self.replace_1)
                            df.loc[mask] = selected_rows
                        else:
                            mask = (df['group_id'] == i - 1) | (df['group_id'] == i)
                            selected_rows = df.loc[mask]
                            selected_rows = selected_rows.iloc[:, :-3].apply(self.replace_1)
                            df.loc[mask] = selected_rows
                    elif i + 1 in group_list:
                        mask = (df['group_id'] == i) | (df['group_id'] == i + 1)
                        selected_rows = df.loc[mask]
                        selected_rows = selected_rows.iloc[:, :-3].apply(self.replace_1)
                        df.loc[mask] = selected_rows
                    elif i - 1 in group_list:
                        mask = (df['group_id'] == i - 1) | (df['group_id'] == i)
                        selected_rows = df.loc[mask]
                        selected_rows = selected_rows.iloc[:, :-3].apply(self.replace_1)
                        df.loc[mask] = selected_rows
                    else:
                        print("合并0采集点时出错，报错group_id = ", i)

                    break

        df['sum'] = df.iloc[:, :-3].sum(axis=1)
        df['group'] = df.apply(self.group_rows, axis=1)
        df["group_id"] = ""
        df.iloc[0, -1] = 1
        n = 1
        m = 2
        while n <= len(df) - 1:
            if df.iloc[n, -2] == df.iloc[n - 1, -2]:
                df.iloc[n, -1] = df.iloc[n - 1, -1]
                n = n + 1
            else:
                df.iloc[n, -1] = m
                n = n + 1
                m = m + 1

        group_list = df['group_id'].unique().tolist()

        # yuannote: 这里要求用户输入SIM分段最大段数
        # sim_sig_max = 99

        while len(group_list) > sim_sig_max:

            result_df = pd.DataFrame(columns=["pattern", "sensitivity_damage", "row_number"])
            i = 1

            while i < df['group_id'].max():
                temp_df = df[(df['group_id'] == i) | (df['group_id'] == i + 1)]
                old_ion_sum = temp_df["sum"].sum()
                temp_df = temp_df.dropna(axis=1, how="all")
                new_ion_sum = temp_df.shape[0] * (temp_df.shape[1] - 3)
                intensity_damage = new_ion_sum - old_ion_sum
                row_num = temp_df.shape[0]
                result_df.loc[len(result_df)] = [i, intensity_damage, row_num]
                i = i + 1

            result_df = result_df.loc[result_df['sensitivity_damage'] ==
                                      result_df['sensitivity_damage'].min()]

            # todo: 现认为此处算法可取得最优结果，但运算过慢
            # 可优化为：每次sensitivity_damage最小的合并方式，如果有多个且sensitivity_damage
            # 数值相等，那么都执行这些合并，但要注意15-16、16-17、17-18的sensitivity_damage都最小
            # 且相等时，只能取15-16、17-18的合并方式，即合并所用的n的list，其中每个n的间隔要大于1
            # 此外运行前要检查，这些合并都执行后，len(group_list)是否小于99，如果小于99，就按照慢方法
            # 一个一个合并
            if result_df.shape[0] > 1:
                result_df.sort_values(by="row_number", inplace=True, ascending=False)
                result_df = result_df[:1]
            else:
                result_df = result_df

            print("*" * 30)
            print(result_df)

            n = result_df.iloc[0, 0]
            temp_df = df[(df['group_id'] == n) | (df['group_id'] == n + 1)]
            temp_df = temp_df.apply(self.replace_1)
            temp_df['sum'] = temp_df.iloc[:, :-3].sum(axis=1)
            temp_df['group'] = temp_df.apply(self.group_rows, axis=1)
            temp_df.loc[:, 'group_id'] = n
            mask = (df['group_id'] == n) | (df['group_id'] == n + 1)
            df.loc[mask, :] = temp_df

            df['group_id'] = df['group_id'].apply(lambda x: x - 1 if x > n else x)

            group_list = df['group_id'].values.tolist()
            group_list = list(set(group_list))
            print(len(group_list))

        df = df[df['sum'] != 0]

        df["dwell_time"] = ""

        # 最小驻留时间，及条件允许情况下每秒最少循环次数在此定义

        # yuannote: 这里要求用户输入最小驻留时间、最少循环次数（if possible）和
        # 实际最少循环次数，默认分别为10、2、0.5，输入0.5后，计算出最大循环时间
        # min_dwell_time = 10

        # min_point_per_s = 2

        # min_point_per_s_limit = 1

        max_cycle_time = (1 / min_point_per_s_limit) * 1000

        for i in df['group_id'].values.tolist():

            ion_num = df.loc[df["group_id"] == i, "sum"].values[0]

            # 一个变量同时满足 if 和 elif 条件时，只有 if 后面的语句会被执行
            # todo: if条件需精简
            if ion_num * min_dwell_time < 1000 / min_point_per_s:
                df.loc[df["group_id"] == i, "dwell_time"] = (1000 / min_point_per_s) / ion_num
                # 如果每秒2个点可满足驻留时间大于最短驻留时间，则取每秒2点下的驻留时间
            elif ion_num * min_dwell_time <= max_cycle_time:
                df.loc[df["group_id"] == i, "dwell_time"] = max_cycle_time / ion_num
                # 如果不满足，查看满足最短驻留时间是否会导致循环时间>1000 ms
                # 若循环时间在1000以内，则优先保驻留时间(灵敏度)
            elif ion_num * 1 <= max_cycle_time:
                df.loc[df["group_id"] == i, "dwell_time"] = max_cycle_time / ion_num
                # 如果不满足，检查循环时间为1000s时，驻留时间是否能大于1，可以的话仍然保驻留时间
            else:
                df.loc[df["group_id"] == i, "dwell_time"] = 1
                # 如果不满足，使其驻留时间=1 ms，此时循环时间上不封顶

        # print(df)
        # yuannote: 这里输出第4个结果
        df.to_excel(outpath + '/SIM_seg_result.xlsx', index=True)

        ##############转安捷伦xml###############
        # yuannote: 这里要求用户选择是否转成安捷伦采集方法，为1则转，默认为1，界面上用打钩形式实现
        # convert_to_ag_method = True

        if convert_to_ag_method:

            from lxml import etree

            group_list = df['group_id'].unique().tolist()
            # 安捷伦方法中的增益因子在此设置
            gain_factor = 10

            mode = "voc"  # mtb or voc

            root = etree.Element('MSAcqMethod', nsmap={'xsi': 'http://www.w3.org/2001/XMLSchema-instance',
                                                       'xsd': 'http://www.w3.org/2001/XMLSchema'})

            msInstrument = etree.SubElement(root, "msInstrument")
            msInstrument.text = ("QQQ")

            ionSource = etree.SubElement(root, "ionSource")
            ionSource.text = ("EI")

            tuneFile = etree.SubElement(root, "tuneFile")
            tuneFile.text = ("atunes.eiex.tune.xml")

            stopMode = etree.SubElement(root, "stopMode")
            stopMode.text = ("ByChromatographTime")

            stopTime = etree.SubElement(root, "stopTime")
            stopTime.text = ("1")

            solventDelay = etree.SubElement(root, "solventDelay")
            solventDelay.text = str(solvent_delay)

            collisionGasOn = etree.SubElement(root, "collisionGasOn")
            collisionGasOn.text = ("true")

            sourceParameters = etree.SubElement(root, "sourceParameters")

            sourceParameter = etree.SubElement(sourceParameters, "sourceParameter")

            para_id = etree.SubElement(sourceParameter, "id")
            para_id.text = ("SourceHeater")

            posPolarityValue = etree.SubElement(sourceParameter, "posPolarityValue")
            posPolarityValue.text = ("250")

            negPolarityValue = etree.SubElement(sourceParameter, "negPolarityValue")
            negPolarityValue.text = ("250")

            isTimeFilterEnabled = etree.SubElement(root, "isTimeFilterEnabled")
            isTimeFilterEnabled.text = ("true")

            timeFilterPeakWidth = etree.SubElement(root, "timeFilterPeakWidth")
            timeFilterPeakWidth.text = ("0.0133333337")

            timeFilter = etree.SubElement(root, "timeFilter")

            activeCount = etree.SubElement(timeFilter, "activeCount")
            activeCount.text = ("1")

            definition = etree.SubElement(timeFilter, "definition")

            time = etree.SubElement(definition, "time")
            time.text = ("0")

            peakWidth = etree.SubElement(definition, "peakWidth")
            peakWidth.text = ("0.0133333337")

            definition_2 = etree.SubElement(timeFilter, "definition")

            time_2 = etree.SubElement(definition_2, "time")
            time_2.text = ("10")

            peakWidth_2 = etree.SubElement(definition_2, "peakWidth")
            peakWidth_2.text = ("0.05")

            useGain = etree.SubElement(root, "useGain")
            useGain.text = ("true")

            if mode == "mtb":

                pass

            else:

                enableNR = etree.SubElement(root, "enableNR")
                enableNR.text = ("true")
                # yuannote: 发现检测代谢样品时，enableNR这行不能出现，反应在方法上为

            timeSegments = etree.SubElement(root, "timeSegments")

            index_n = 0
            for i in group_list:

                timeSegment = etree.SubElement(timeSegments, "timeSegment")

                index = etree.SubElement(timeSegment, "index")
                index_n += 1
                index.text = str(int(index_n))

                startTime = etree.SubElement(timeSegment, "startTime")
                startTime.text = "{:.4f}".format(df.index[df["group_id"] == i][0])

                isDataSaved = etree.SubElement(timeSegment, "isDataSaved")
                isDataSaved.text = ("true")

                scanSegments = etree.SubElement(timeSegment, "scanSegments")

                scanSegment = etree.SubElement(scanSegments, "scanSegment")

                seg_index = etree.SubElement(scanSegment, "index")
                seg_index.text = ("1")

                scanType = etree.SubElement(scanSegment, "scanType")
                scanType.text = ("MS2SIM")

                scanElements = etree.SubElement(scanSegment, "scanElements")

                temp_df = df[(df['group_id'] == i)].copy()
                dwell_time = round(temp_df.iloc[0, -1], 1)
                temp_df = temp_df.dropna(axis=1, how="all")
                ion_list = list(temp_df.iloc[:, :-4].columns)
                ion_list.sort(reverse=True)

                for n, ion in enumerate(ion_list, start=1):
                    scanElement = etree.SubElement(scanElements, "scanElement")

                    ele_index = etree.SubElement(scanElement, "index")
                    ele_index.text = str(n)

                    compoundName = etree.SubElement(scanElement, "compoundName")
                    compoundName.text = ("seg_{}_ion_{}".format(i, ion))

                    isISTD = etree.SubElement(scanElement, "isISTD")
                    isISTD.text = ("false")

                    ms2LowMz = etree.SubElement(scanElement, "ms2LowMz")
                    ms2LowMz.text = (str(ion))

                    ms2Res = etree.SubElement(scanElement, "ms2Res")
                    ms2Res.text = ("LowRes")

                    dwell = etree.SubElement(scanElement, "dwell")
                    dwell.text = str(dwell_time)

                    gain = etree.SubElement(scanElement, "gain")
                    gain.text = str(gain_factor)

            instrumentCurves = etree.SubElement(root, "instrumentCurves")

            samplingRate = etree.SubElement(instrumentCurves, "samplingRate")
            samplingRate.text = ("5")

            chromatograms = etree.SubElement(root, "chromatograms")

            chromatogram = etree.SubElement(chromatograms, "chromatogram")

            chrom_index = etree.SubElement(chromatogram, "index")
            chrom_index.text = ("1")

            chromType = etree.SubElement(chromatogram, "chromType")
            chromType.text = ("TIC")

            label = etree.SubElement(chromatogram, "label")
            label.text = ("TIC")

            offset = etree.SubElement(chromatogram, "offset")
            offset.text = ("0")

            yRange = etree.SubElement(chromatogram, "yRange")
            yRange.text = ("1E+07")

            defaultExtractionWindow = etree.SubElement(chromatograms, "defaultExtractionWindow")

            minus = etree.SubElement(defaultExtractionWindow, "minus")
            minus.text = ("0.3")

            plus = etree.SubElement(defaultExtractionWindow, "plus")
            plus.text = ("0.7")

            etree.indent(root, space="  ", level=0)

            # yuannote: 如果要转安捷伦方法，这里输出第5个结果
            with open(outpath + "/qqqacqmethod.xml", "wb") as f:
                f.write(etree.tostring(root, encoding="UTF-8", xml_declaration=True))

        print("all done")

# msp_path = r"C:\Users\86724\Desktop\采集方法测试结果\Remove_Duplicates.msp"
# rt_data_path = r"C:\Users\86724\Desktop\采集方法测试结果\New_RT_list.xlsx"
# set_name_list = True
# name_list_path = r"C:\Users\86724\Desktop\采集方法测试结果\toamto.txt"
# mz_min = 35
# mz_max = 400
# outpath = r"C:\Users\86724\Desktop"
# rt_window = 0.5
# min_ion_intensity_percent = 5
# min_ion_num = 2
# prefer_mz_threshold = 60
# similarity_threshold = 0.85
# fr_factor = 2
# retention_time_max = 68.8
# solvent_delay = 0
# sim_sig_max = 99
# min_dwell_time = 10
# min_point_per_s = 2
# min_point_per_s_limit = 1
# convert_to_ag_method = True
#
#
# a = GetMethod()
# a.Main(msp_path, rt_data_path, set_name_list, name_list_path, mz_min, mz_max, outpath, rt_window,
#              min_ion_intensity_percent, min_ion_num, prefer_mz_threshold, similarity_threshold, fr_factor, retention_time_max, solvent_delay, sim_sig_max, min_dwell_time,
#              min_point_per_s, min_point_per_s_limit, convert_to_ag_method)
# print("all done")
