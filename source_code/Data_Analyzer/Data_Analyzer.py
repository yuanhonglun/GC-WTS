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
import netCDF4 as nc

class DataAnalysis():

    def find_nearest(self, array, value):
        """
        Find the index of the nearest element to a given value in a sorted array.

        Args:
            array (numpy.ndarray): The sorted array in which to find the nearest element.
            value (float): The value to which we want to find the nearest element.

        Returns:
            int: The index of the nearest element in the array.
        """
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or abs(value - array[idx - 1]) < abs(value - array[idx])):
            return idx - 1
        else:
            return idx

    def calulate_lwma(self, x, weights, no_smooth_number):
        """
        Calculate the weighted moving average (LWMA) for a given series 'x' using provided 'weights'.

        Args:
            x (pandas.Series): The series for which to calculate LWMA.
            weights (list): The list of weights for the LWMA calculation.
            no_smooth_number (int): The number of data points for smoothing.

        Returns:
            float: The calculated LWMA value.
        """

        if isinstance(x.index[0], str):
            list_sum = sum(weights[no_smooth_number - int(x.index[0]):])
            return np.dot(x, weights) / list_sum
        elif isinstance(x.index[-1], str):
            list_sum = sum(weights[int(x.index[-1]) + 1:])
            return np.dot(x, weights) / list_sum
        else:
            list_sum = sum(weights)
            return np.dot(x, weights) / list_sum

    def validateTitle(self, title):
        """
        Validate and sanitize a title string by replacing special characters with underscores.

        Args:
            title (str): The input title string to be validated.

        Returns:
            str: The sanitized title string with special characters replaced by underscores.
        """
        rstr = r"[\/\\\:\*\?\"\<\>\|]"  # '/ \ : * ? " < > |'
        new_title = re.sub(rstr, "_", title)  # 替换为下划线
        return new_title

    def lwma(self, local_df, n):
        """
        Calculate the Linear Weighted Moving Average (LWMA) for each column in a DataFrame.

        Args:
            local_df (pd.DataFrame): The input DataFrame.
            n (float): The window size for LWMA calculation.

        Returns:
            pd.DataFrame: The DataFrame containing LWMA values for each column.
        """
        n = int(n)
        smooth_df = pd.DataFrame(columns=local_df.columns.values, index=local_df.index.values)
        list1 = [i for i in range(1, int((n + 1) / 2) + 1)]
        weights = list1 + list(reversed(list1))[1:]
        no_smooth_number = int((n - 1) / 2)
        tmp_df = pd.DataFrame(np.zeros((no_smooth_number, local_df.shape[1])), columns=list(local_df.columns))
        tmp_df = tmp_df.set_axis(list(map(str, list(tmp_df.index))))
        local_df_padding = pd.concat([tmp_df, local_df, tmp_df])
        for ion in local_df.columns.values:
            smooth = local_df_padding[ion].rolling(window=n, min_periods=n, center=True).apply(self.calulate_lwma,
                                                                                               args=(
                                                                                                   weights,
                                                                                                   no_smooth_number)).to_list()
            smooth_df[ion] = smooth[no_smooth_number:-no_smooth_number]

        smooth_df = smooth_df.dropna(axis=0, how='all')

        return smooth_df

    def derivative(self, smooth_df, m):
        """
        Calculate derivatives and noise characteristics for each column in a DataFrame.

        Args:
            smooth_df (pd.DataFrame): The input DataFrame containing smoothed data.
            m (float): A threshold value for noise filtering.

        Returns:
            smooth_df (pd.DataFrame): The DataFrame containing derivatives and noise characteristics.
            noise_df (pd.DataFrame): The DataFrame containing noise characteristics.
        """
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

            ff = np.median(np.asarray([i for i in abs_fd_1 if 0 < i < max(abs_fd_1) * 0.1]))
            sf = np.median(np.asarray([i for i in neg_sd_1 if i < 0 < abs(i) and abs(i) < max(neg_sd_1) * 0.1]))
            af = np.median(np.asarray([i for i in abs_ad_1 if 0 < i < max(abs_ad_1) * 0.1]))

            ff = max(ff, m)
            if af < m:
                af = m

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

    def find_peak(self, ion_list, smooth_df, noise_df, filter_factor, df):
        """
        Find peaks and perform peak-related calculations.

        Args:
            ion_list (list): List of ions to process.
            smooth_df (pd.DataFrame): The input DataFrame containing smoothed data.
            noise_df (pd.DataFrame): The DataFrame containing noise characteristics.
            filter_factor (float): A filtering factor for peak detection.

        Returns:
            peak_dic (dict): A dictionary containing peak information for each ion.
            con_peak_dic (dict): A dictionary containing continuous peaks for each ion.
            smooth_df_final (pd.DataFrame): The final smoothed DataFrame with baseline correction.
        """
        peak_dic = {}
        con_peak_dic = {}
        NF = self.pre_de_redundancy_peak(smooth_df, ion_list)
        for ion in ion_list:
            peak_df = pd.DataFrame(columns=["left", "apex", "right", "SV"])
            ff = noise_df.loc["ff", ion]
            sf = noise_df.loc["sf", ion]
            f_list = smooth_df[str(ion) + "_fd"].values.tolist()
            s_list = smooth_df[str(ion) + "_sd"].values.tolist()
            a_list = smooth_df[str(ion) + "_ad"].values.tolist()
            rt_list = smooth_df.index.values
            i_list = smooth_df[ion].values.tolist()
            m = 0
            total_list = []
            index_list = []
            aa = 1
            while m <= len(f_list) - 3:
                peak_list = []
                if aa != m:
                    aa = m
                    if f_list[m] > ff * filter_factor and f_list[m + 1] > ff * filter_factor and i_list[m] > 0 and \
                            i_list[
                                m + 1] > 0:
                        w = m - 2
                        if w < 0:
                            w = 0
                        tmp_list = i_list[w: m + 3]
                        u = tmp_list.index(i_list[m])
                        q = tmp_list.index(min(tmp_list))
                        p = m - (u - q)
                        if i_list[p] == 0:
                            p = m
                            left = rt_list[p]
                        else:
                            left = rt_list[p]

                        peak_list.append(left)

                        for x in range(p + 1, len(f_list) - 3):
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
                                        break

                                    elif i_list[y] < i_list[x] * 0.05:
                                        if rt_list[y] <= rt_list[x]:
                                            peak_list.append(rt_list[y + 1])
                                            m = y + 1
                                        else:
                                            peak_list.append(rt_list[y])
                                            m = y
                                        break
                                break
                            elif x > len(f_list) - 2:

                                peak_list = []
                                m = x
                                break

                            else:
                                if f_list[x - 1] > 0 and s_list[x] < -1 * sf:
                                    if f_list[x] < 0 or f_list[x + 1] < 0:
                                        start_index = max(0, x - 2)
                                        end_index = min(len(i_list), x + 3)
                                        five_elements = i_list[start_index:end_index]
                                        max_index = five_elements.index(max(five_elements))
                                        n = start_index + max_index
                                        peak_list.append(rt_list[n])
                                        t = smooth_df.loc[rt_list[n], ion]
                                        index_list.append("peak_" + str(rt_list[n]))

                                        for y in range(n + 1, len(f_list) - 3):
                                            if f_list[y] > -ff * filter_factor and f_list[y + 1] > -ff * filter_factor:
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
                                                elif rt_list[p_1] <= rt_list[n]:
                                                    right = rt_list[y]
                                                    m = y
                                                else:
                                                    right = rt_list[p_1]
                                                    m = p_1
                                                peak_list.append(right)
                                                break

                                            elif i_list[y] < i_list[n] * 0.05:
                                                if rt_list[y] <= rt_list[n]:
                                                    peak_list.append(rt_list[y + 1])
                                                    m = y + 1
                                                else:
                                                    peak_list.append(rt_list[y])
                                                    m = y
                                                break
                                            elif y == len(f_list) - 4:
                                                peak_list.append(rt_list[y])
                                                m = y
                                                break

                                        break
                        if peak_list != [] and len(peak_list) == 3:
                            total_list.append(peak_list)
                        elif len(peak_list) == 2:
                            index_list.remove("peak_" + str(rt_list[n]))
                    else:
                        m += 1
                else:
                    break
            peak_df = pd.DataFrame(total_list, columns=["left", "apex", "right"], index=index_list)
            peak_df = self.trailing_peak_filtration(peak_df, df, ion)
            fb = 'smooth_df_final' in dir()
            if fb == False:
                from copy import copy
                smooth_df_final = copy(smooth_df.applymap(lambda x: x if x >= 0 else 0))
            smooth_df_final = self.baseline_correction(peak_df, smooth_df, ion, smooth_df_final)
            t1 = list(smooth_df_final.loc[:, ion])
            t2 = list(smooth_df_final.index)
            peak_df = self.de_redundancy_peak(smooth_df_final, smooth_df, peak_df, NF, ion)
            rt_list = rt_list.tolist()
            for peak in peak_df.index.values:
                apex_rt = float(peak_df.loc[peak, "apex"])
                apex_i = float(smooth_df_final.loc[apex_rt, ion])
                apex_rt_index = rt_list.index(apex_rt)

                apex_left_rt = float(rt_list[apex_rt_index - 1])
                apex_left_i = float(smooth_df_final.loc[apex_left_rt, ion])
                apex_right_rt = float(rt_list[apex_rt_index + 1])
                apex_right_i = float(smooth_df_final.loc[apex_right_rt, ion])

                new_apex_rt, new_apex_i = self.second_order_fit_X(
                    apex_left_rt, apex_left_i, apex_rt, apex_i, apex_right_rt, apex_right_i)
                smooth_df_final.loc[apex_rt, ion] = np.nan
                smooth_df_final.loc[new_apex_rt, ion] = new_apex_i
                smooth_df_final.sort_index(inplace=True)
                peak_df.loc[peak, "apex"] = new_apex_rt
                left_rt = peak_df.loc[peak, "left"]
                left_i = smooth_df_final.loc[left_rt, ion]
                right_rt = peak_df.loc[peak, "right"]
                right_i = smooth_df_final.loc[right_rt, ion]
                x_df = smooth_df_final.loc[left_rt: right_rt, ion]
                x_df = x_df.dropna(axis=0, how="any")
                y = x_df.values
                del x_df, apex_left_rt, apex_right_rt
                gc.collect()
            temp_df = smooth_df_final[[ion]]
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
            con_peak_dic[ion] = con_peak_df
            peak_dic[ion] = peak_df
        return peak_dic, con_peak_dic, smooth_df_final

    def trailing_peak_filtration(self, peak_df, raw_df, ion):
        """
        Filter trailing peaks and identify incorrectly located peaks.

        Args:
            peak_df (pd.DataFrame): DataFrame containing peak information.
            raw_df (pd.DataFrame): The input DataFrame containing raw data.
            ion (str): Ion identifier.

        Returns:
            peak_df (pd.DataFrame): DataFrame with filtered peaks.
        """
        for peak in peak_df.index.values:
            try:
                apex_rt = peak_df.loc[peak, "apex"]
                apex_i = raw_df.loc[apex_rt, ion]
                left_rt = peak_df.loc[peak, "left"]
                left_i = raw_df.loc[left_rt, ion]
                left_left_i = raw_df.shift(1).loc[left_rt, ion]
                right_left_i = raw_df.shift(-1).loc[left_rt, ion]
                right_rt = peak_df.loc[peak, "right"]
                right_i = raw_df.loc[right_rt, ion]
                right_right_i = raw_df.shift(-1).loc[right_rt, ion]
                left_right_i = raw_df.shift(1).loc[right_rt, ion]
                if left_rt < apex_rt < right_rt:
                    if right_i > 0.85 * apex_i or left_i > 0.85 * apex_i:
                        peak_df = peak_df.drop(peak, axis=0)
                    elif right_left_i == 0 or left_left_i == 0 or left_right_i == 0 or right_right_i == 0:
                        peak_df = peak_df.drop(peak, axis=0)
            except KeyError:
                pass

        return peak_df

    def second_order_fit_X(self, x1, y1, x2, y2, x3, y3):
        if x1 < x2 < x3 and y2 > y1 and y2 > y3:

            X = np.array([x1, x2, x3])
            Y = np.array([y1, y2, y3])
            coef = np.polyfit(X, Y, 2)
            x2_new = (-coef[1]) / (2 * coef[0])

            if x1 < x2_new < x3:
                y2_new = coef[0] * (x2_new ** 2) + coef[1] * x2_new + coef[2]
                if y2_new > y1 and y2_new > y3 and y2_new <= 1.5 * y2:
                    y2_new = y2_new
                else:
                    y2_new = y2
            else:
                x2_new = x2
                y2_new = y2
            return x2_new, y2_new

        else:

            return x2, y2

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
        """
        Deconvolve peaks in the ion list.

        Args:
            ion_list (list): List of ions to process.
            smooth_df_final (pd.DataFrame): Smoothed data.
            con_peak_dic (dict): Dictionary of con_peak DataFrames.

        Returns:
            dict: Dictionary containing deconvolved peaks for each ion.
            pd.DataFrame: Deconvolved data.
        """
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
                    para_list = []
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
                    before_half_y = y_0[:(len(y_0) // 2)]
                    after_half_y = y_0[(len(y_0) // 2):]
                    background_min = max(min(y_0), 0)
                    background_max = max(min(before_half_y), min(after_half_y), 0.01)
                    if background_min < background_max:
                        para_list.append(background_min)
                        bound_min.append(background_min)
                        bound_max.append(background_max)
                    else:
                        para_list.append(0)
                        bound_min.append(0)
                        bound_max.append(1000000000)

                    popt, pcov = curve_fit(self.Gaussian, x, y, bounds=(bound_min, bound_max), p0=para_list,
                                           maxfev=100000)
                    fit = self.Gaussian(x, *popt)
                    y_list = self.fit_plot(x, *popt)
                    baseline = np.zeros_like(x) + popt[-1]
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
                                df = pd.DataFrame(peak, index=x, columns=[peak_name])
                                decon_data_df = decon_data_df.combine_first(df)
        return decon_peak_dic, decon_data_df


    def sharpness(self, peak_df, data_df):
        data_df["number_index"] = [x for x in range(len(data_df))]
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
                sharpness = max(left_list)
            elif right_list:
                sharpness = max(right_list)
            else:
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
            peak_dic[ion] = self.sharpness(peak_dic[ion], data_df)

        for ion in decon_peak_dic.keys():
            new_peak_df = pd.DataFrame(columns=["left", "apex", "right", "SV"])
            for peak in decon_peak_dic[ion].index.values:
                data_df = decon_data_df.loc[:, [str(peak)]]
                data_df = data_df.dropna(axis=0, how="any")
                peak_df = decon_peak_dic[ion].loc[[peak], :]
                peak_df = self.sharpness(peak_df, data_df)
                new_peak_df = pd.concat((new_peak_df, peak_df))
            decon_peak_dic[ion] = new_peak_df
        return peak_dic, decon_peak_dic

    def group_peak(self, peak_dic, decon_peak_dic, smooth_df, smooth_df_final, decon_data_df, wid, sigma, bin_num,
                   group_method):
        """
        Group peaks based on specified criteria.

        Args:
            peak_dic (dict): Dictionary containing peak information.
            decon_peak_dic (dict): Dictionary containing deconvolved peak information.
            smooth_df (pd.DataFrame): Smoothed data.
            smooth_df_final (pd.DataFrame): Final smoothed data.
            decon_data_df (pd.DataFrame): Deconvolved data.
            wid (float): Width parameter.
            sigma (float): Sigma parameter.
            bin_num (int): Number of bins.
            group_method (int): Grouping method identifier (0 for MSDIAL, 1 for AMDIS).

        Returns:
            pd.DataFrame: Matched wave data.
            pd.DataFrame: Peak group data.
            pd.DataFrame: Quantitative results.
        """
        quant_result_df = pd.DataFrame(columns=['Peak_group', 'Ion', "Relative_Peak_Area", "Peak_Height"])
        rt_list = smooth_df.index.values.tolist()
        if bin_num >= 1:
            bin_num = round(bin_num)
            new_list = [rt_list[0]] + [rt_list[i] + (rt_list[i + 1] - rt_list[i]) / bin_num * j for i in
                                       range(len(rt_list) - 1)
                                       for j in range(1, bin_num + 1)] + [rt_list[-1]]
        elif 0 < bin_num < 1:
            bin_num_re = math.ceil(1 / bin_num)
            new_list = [rt_list[0]]
            idx = bin_num_re
            while idx < len(rt_list) - 1:
                new_list.append(rt_list[idx])
                idx += bin_num_re
            new_list.append(rt_list[-1])
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
                        if matched_wave.iloc[i, 2] == 0:
                            matched_wave.iloc[i, 2] = left
                        elif matched_wave.iloc[i, 2] > left:
                            matched_wave.iloc[i, 2] = left
                        if matched_wave.iloc[i, 3] == 0:
                            matched_wave.iloc[i, 3] = right
                        elif matched_wave.iloc[i, 3] < right:
                            matched_wave.iloc[i, 3] = right
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
        if group_method == 0:
            matched_wave, peak_group_df = self.MSDIAL(matched_wave, wid, sigma, peak_dic, decon_peak_dic)
        elif group_method == 1:
            matched_wave, peak_group_df = self.AMDIS(matched_wave)
        for group, row in peak_group_df.iterrows():
            quant_df = row["ions"]
            quant_df_index = quant_df.index.tolist()
            Quant_Ion = int(quant_df_index[0])
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
                    sliced_df = decon_data_df.loc[
                                (decon_data_df.index >= row["left"]) & (decon_data_df.index <= row["right"]),
                                :]
                    sliced_df = sliced_df.filter(like=f"{Quant_Ion}_peak_{col}")
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
        return matched_wave, peak_group_df, quant_result_df

    def peak_group_df_add_retention_infor(self, peak_group_df, RI_Presentation, standard_df, RI_min, RI_max):
        """
        peak_group_df: Peak group information
        RI_Presentation: RI calculation mode
        standard_df: Reference standard table for RI calculation
        RI_min, RI_max: RI thresholds
        """
        for i, row in peak_group_df.iterrows():
            rt_sample = row['rt'] / 60
            if RI_Presentation == "Kovats":
                peak_group_df.loc[i, "Kovats_RI"] = self.RT_to_Kovats_RI_transform(rt_sample, standard_df, RI_min,
                                                                                   RI_max)
            elif RI_Presentation == "Fiehn":
                peak_group_df.loc[i, "Fiehn_RI"] = self.RT_to_Fiehn_RI_transform(rt_sample, standard_df)
        return peak_group_df

    def MSDIAL(self, matched_wave, wid, sigma, peak_dic, decon_peak_dic):

        sg_list = []
        for n in range(0, wid * 2 + 1):
            sg_list.append((1 - (pow(((-wid + n) / sigma), 2))) * (math.exp(-0.5 * pow(((-wid + n) / sigma), 2))))
        matched_wave["gauss_SV"] = matched_wave["SV"].rolling(window=wid * 2 + 1, min_periods=wid * 2 + 1,
                                                              center=True).apply(
            lambda x: sum(np.multiply(np.asarray(x), np.asarray(sg_list)))).to_list()
        matched_wave = matched_wave.dropna(axis=0, how="any")
        peak_group_df = pd.DataFrame(columns=matched_wave.columns)
        for n in range(1, matched_wave.shape[0] - 2):  ##避免出现sv顶点在最后一个情况
            if matched_wave.iloc[n, 5] > 0 and matched_wave.iloc[n, 5] - matched_wave.iloc[n - 1, 5] > 0 and \
                    matched_wave.iloc[n, 5] - matched_wave.iloc[n + 1, 5] > 0:
                n_left = n
                n_right = n
                while matched_wave.iloc[n_right, 5] >= matched_wave.iloc[n_right + 1, 5]:
                    n_right += 1
                    ##避免右点刚好是最后一个点，n_right + 1索引错误
                    if n_right >= matched_wave.shape[0] - 1:
                        break
                while matched_wave.iloc[n_left, 5] >= matched_wave.iloc[n_left - 1, 5]:
                    n_left -= 1
                tmp_ions_df = pd.concat(matched_wave.iloc[n_left:n_right + 1, [4]]["ions"].tolist())
                tmp_ions_df.sort_values('intensity', ascending=False, inplace=True)
                duplicate_index = tmp_ions_df.index.duplicated()
                tmp_ions_df = tmp_ions_df[~duplicate_index]
                matched_wave.at[matched_wave.index[n], 'ions'] = tmp_ions_df
                peak_group_df.loc[len(peak_group_df)] = matched_wave.iloc[n]

        peak_group_df = peak_group_df.drop(peak_group_df[peak_group_df['SV'] == 0].index)

        return matched_wave, peak_group_df

    def AMDIS(self, matched_wave):

        matched_wave.insert(5, "compound_point", np.nan, allow_duplicates=False)
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
        return matched_wave, peak_group_df

    def is_number(self, s):

        pattern = r'^[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?$'
        return bool(re.match(pattern, s))

    def remove_content_from_list(self, lst, content):
        pattern = re.escape(content)  # 转义特殊字符
        regex = re.compile(pattern)
        return [regex.sub('', string) for string in lst if isinstance(string, str)]

    def read_msp(self, msp):
        """
        Read data from an MSP file and parse it into a dictionary.

        Args:
            msp (str): The path to the MSP file.

        Returns:
            dict: A compounds' dictionary containing parsed data from the MSP file.
            pd.DataFrame: Compounds' RI information.
        """
        msp_file = open(msp, "r")
        list_1 = msp_file.readlines()
        new_list = [item.replace('NAME:', 'Name:') for item in list_1]
        list_1 = [item.replace("Num peaks:", "Num Peaks:") for item in new_list]
        lists = str(list_1)
        lines = lists.split("Name: ")
        meta = {}
        for l in lines:
            line1 = l.strip().split("\\n")
            name_1 = line1[0]
            line2 = l.strip().split("Num Peaks:")
            ion_intens_dic = {}
            if len(line2) > 1:
                if ';' in line2[1]:
                    matches = re.findall(r"(\d+) (\d+);", line2[1])
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
                        if len(ion1) == 2 and self.is_number(ion1[0]) and self.is_number(ion1[1]):
                            key = round(float(ion1[0]))
                            value = float(ion1[1])
                            if key in ion_intens_dic:
                                ion_intens_dic[key] = max(ion_intens_dic[key], value)
                            else:
                                ion_intens_dic[key] = value
                elif '\\n' in line2[1]:
                    line3 = line2[1].split("\\n', '")[1:-2]
                    ion_intens_dic = {}
                    for ion in line3:
                        ion1 = ion.split(" ")
                        if len(ion1) == 2 and self.is_number(ion1[0]) and self.is_number(ion1[1]):
                            key = round(float(ion1[0]))
                            value = float(ion1[1])
                            if key in ion_intens_dic:
                                ion_intens_dic[key] = max(ion_intens_dic[key], value)
                            else:
                                ion_intens_dic[key] = value
                        elif ':' in ion:
                            pattern = re.compile(r'(\d+):(\d+)')
                            matches = pattern.findall(ion)
                            for key, value in matches:
                                key = round(float(key))
                                value = int(value)
                                if key in ion_intens_dic:
                                    ion_intens_dic[key] = max(ion_intens_dic[key], value)
                                else:
                                    ion_intens_dic[key] = value

                else:
                    print('The format is not recognized.')
            meta[name_1] = ion_intens_dic
        RI_df = self.read_msp_RI(msp)
        return meta, RI_df

    def group_cmp_inf(self, lines):
        """
        FUNCTION: Retrieve the line containing the substance name and translate it into English in academic paper style
        """
        group_inf_idx = [i for i, p in enumerate(lines) if 'Name:' in p]
        group_inf_idx.append(len(lines))
        return group_inf_idx

    def read_msp_RI(self, msp):
        """
        FUNCTION: Convert Sample RT-RI
        RT_data: RT library
        msp: MSP library
        return: msp_RI
        """
        msp_file = open(msp, "r")
        lines = msp_file.readlines()
        group_inf_idx = self.group_cmp_inf(lines)
        RI_df = pd.DataFrame(columns=['Name', 'RI'])
        for j in range(len(group_inf_idx) - 1):
            group_inf = lines[group_inf_idx[j]:group_inf_idx[j + 1]]
            prefixes = [r'SemiStdNP=\d+', r'RI:\d+\n', r'Any=\d+', r'RETENTIONINDEX: \d+\.\d+', r'Synon: RI: \d+']
            pattern = "|".join(prefixes)
            for string in group_inf:
                if 'Name:' in string:
                    RI_list = [string.replace('Name: ', '')]
                if matches := re.findall(pattern, string):
                    RI_list.extend([int(re.findall(r"\d+", match)[0]) for match in matches])
                    RI_df.loc[len(RI_df.index)] = RI_list
        RI_df['Name'] = RI_df['Name'].str.rstrip('\n')
        RI_df.set_index('Name', inplace=True)
        return RI_df

    def dot_product_distance(self, p, q):
        if (np.sum(p)) == 0 or (np.sum(q)) == 0:
            score = 0
        else:
            score = np.power(np.sum(q * p), 2) / \
                    (np.sum(np.power(q, 2)) * np.sum(np.power(p, 2)))
        return score

    def weighted_dot_product_distance(self, compare_df):
        m_q = pd.Series(compare_df.index)
        i_q = np.array(compare_df.iloc[:, 0])
        i_r = np.array(compare_df.iloc[:, 1])
        k = 0.5
        l = 2
        w_q = np.power(i_q, k) * np.power(m_q, l)
        w_r = np.power(i_r, k) * np.power(m_q, l)
        ss = self.dot_product_distance(w_q, w_r)
        shared_spec = np.vstack((i_q, i_r))
        shared_spec = pd.DataFrame(shared_spec)
        shared_spec = shared_spec.loc[:, (shared_spec != 0).all(axis=0)]
        m = int(shared_spec.shape[1])
        if m >= 2:
            FR = 0
            for i in range(1, m):
                s = (shared_spec.iat[0, i] / shared_spec.iat[0, (i - 1)]) * (
                        shared_spec.iat[1, (i - 1)] / shared_spec.iat[1, i])
                if s > 1:
                    s = 1 / s
                FR = FR + s
            ave_FR = FR / (m - 1)
            NU = int(len(compare_df))
            composite_score = ((NU * ss) + (m * ave_FR)) / (NU + m)
        else:
            composite_score = ss

        return composite_score

    def euclidean_intensity(self, x, y, k, b):
        """y = kx + b """
        return math.fabs(b + k * (x - 0) - y) / (math.sqrt(k * k + 1))

    def line_between_dots(self, x1, y1, x2, y2):
        """y = kx + b 返回 k,b"""
        if x1 == x2:
            return f'input_error: x1 = x2 = {x1}', np.nan
        else:
            return (y1 - y2) / (x1 - x2), (x1 * y2 - x2 * y1) / (x1 - x2)

    def line(self, x, k, b):
        return k * x + b

    def base_line(self, rt_list, intensity_list):
        para, pcov = curve_fit(self.line, rt_list, intensity_list)
        y_fitted = self.line(rt_list, para[0], para[1])
        return para[0], para[1], y_fitted

    def baseline_correction(self, peak_df, smooth_df, ion, smooth_df_final):
        base_line_rt = []
        base_line_intensity = []
        list_now = []
        for peak_line in peak_df.itertuples():
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
            son_df["euclidean_intensity_list"] = euclidean_intensity_list
            list_now.extend(list(son_df[ion]))
            son_son_df = son_df.nsmallest(int(son_df.shape[0] / 2), "euclidean_intensity_list", keep='all')
            k1, b1, y_fitted = self.base_line(son_son_df.index, son_son_df.iloc[:, 0])
            for rt1, intensity1 in son_df.iterrows():
                new_intensity = intensity1.values[0] - self.line(rt1, k1, b1)
                if new_intensity >= 0:
                    base_line_intensity.append(new_intensity)
                else:
                    base_line_intensity.append(0)
                base_line_rt.append(rt1)
        for index1, line1 in smooth_df.iterrows():
            if index1 in base_line_rt:
                smooth_df_final.loc[index1, ion] = base_line_intensity[base_line_rt.index(index1)]
        return smooth_df_final

    def pre_de_redundancy_peak(self, smooth_df, ion_list):
        """
        Remove redundant peaks
        1. First extract the [ion] dataframe from smooth_df, then iterate through the 13 peaks.
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
                    if number >= 7:
                        median_value = median(window_list)
                        NF_list.append(abs(median_value - window_mean) / (window_mean ** 0.5))
                m = m + 13
            all_NF_list.extend(NF_list)
        NF = median(all_NF_list)

        return NF

    def de_redundancy_peak(self, smooth_df_final, smooth_df, peak_df, NF, ion):
        index_list = []
        for index, row in peak_df.iterrows():
            tmp_apex_intensity = smooth_df.loc[row["apex"], ion]
            ex_tmp_apex_intensity = smooth_df_final.loc[row["apex"], ion]
            if (4 * NF * ((tmp_apex_intensity) ** 0.5)) > ex_tmp_apex_intensity:
                index_list.append(index)
            elif ex_tmp_apex_intensity == 0:
                #为了去除基线矫正后顶点丰度刚好为0的峰
                index_list.append(index)
        peak_df.drop(index_list, axis=0, inplace=True)

        return peak_df

    ############################

    ############################
    def RT_similarity(self, rt_s, rt_l, window, level_factor, max_penalty):
        """
        The method used in NIST for comparing RI similarity
        rt_s: Sample RT
        rt_l: RT in the library
        penalty_rate: Penalty rate. weak-20, average-50, strong-100, very strong-200
        rt_tolerance: RT tolerance
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
        """
        Suggestion: Possible new compound
        row: Row information traversed in peak_group_df
        score_df: Dataframe of matching scores for compounds in the library
        min_ions_number: Minimum number of ions contained in the signal
        return: results_df annotated with new compounds
        """
        if len(row["ions"]) >= min_ions_number:
            if max(score_df["score"]) < limit_of_score:
                results_df.loc[group_rt_s, "New_Compound"] = ("This compound might be a new compound not present in current library")
        return results_df

    def species_pool_opt(self, results_df, species_pool, species_pool_coefficient=0.3,
                         all_match_optimization_coefficient=0.2):
        """
        Optimize qualitative results based on input species library information.
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
        FUNCTION: Generate mass spectrometry comparison plot
        group_rt_s: Group names
        compound: Compound name in the library
        compare_df_1: Dataframe containing two mass spectrometry information
        """

        compare_df_new = compare_df_1[compare_df_1 >= 10]
        compare_df_new.dropna(inplace=True)

    def match_compare_df_correct(self, match_compare_df):
        intensity_s_sum = match_compare_df["intensity_s"].sum()
        w = 1 / (-1 + 0.5 + intensity_s_sum / 999)
        for index, row in match_compare_df.iterrows():
            A = row[0] / 999
            match_compare_df.loc[index, 'intensity_s'] = A / (1 + w * A)
        return match_compare_df

    def score_correct(self, compare_df_1, r_match_compare_df, MS_score, R_MS_score):

        if compare_df_1.shape[0] == 1:
            MS_score = MS_score * 0.6
        elif compare_df_1.shape[0] == 2:
            MS_score = MS_score * 0.7
        elif compare_df_1.shape[0] == 3:
            MS_score = MS_score * 0.94
        elif compare_df_1.shape[0] == 4:
            MS_score = MS_score * 0.97
        else:
            MS_score = MS_score

        if r_match_compare_df.shape[0] == 1:
            R_MS_score = R_MS_score * 0.75
        elif r_match_compare_df.shape[0] == 2:
            R_MS_score = R_MS_score * 0.88
        elif r_match_compare_df.shape[0] == 3:
            R_MS_score = R_MS_score * 0.94
        elif r_match_compare_df.shape[0] == 4:
            R_MS_score = R_MS_score * 0.97
        else:
            R_MS_score = R_MS_score

        return MS_score, R_MS_score

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

        """
        Perform peak group identification based on provided parameters.

        Args:
            peak_group_df (pd.DataFrame): DataFrame containing peak group information.
            smooth_df (pd.DataFrame): DataFrame with smoothed data.
            smooth_df_final (pd.DataFrame): Final smoothed data.
            meta (dict): Metadata containing ion intensity information.
            RT_df (pd.DataFrame): DataFrame with retention time data.
            decon_peak_dic: Deconvoluted peak data.
            decon_data_df: Deconvoluted data DataFrame.
            peak_group_search: Parameter for peak group search.
            match_weight: Weight for matching.
            r_match_weight: Weight for reverse matching.
            search_wid: Search width.
            retention_score_mode: Retention score mode.
            window: Window parameter.
            window_scale: Scaling factor for the window.
            level_factor: Level factor for scoring.
            max_penalty: Maximum penalty for scoring.
            no_info_penalty: Penalty for missing information.
            minimum_number_of_ions: Minimum number of ions for scoring.
            sim_threshold: Similarity threshold.
            group_weight: Weight for group matching.
            direct_weight: Weight for direct matching.
            inaccurate_ri_threshold: Threshold for inaccurate retention index (RI).
            inaccurate_ri_level_factor: Level factor for inaccurate RI.

        Returns:
            pd.DataFrame: DataFrame containing the results of peak group identification.
        """

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
            if retention_score_mode == "None":
                compound_list = RT_df.index.values.tolist()
            elif retention_score_mode == "RI":
                if 'Kovats_RI' in peak_group_df.columns:
                    if not isinstance(peak_group_df.loc[group, "Kovats_RI"], str):
                        group_ri = float(peak_group_df.loc[group, "Kovats_RI"])
                elif 'Fiehn_RI' in peak_group_df.columns:
                    if not isinstance(peak_group_df.loc[group, "Fiehn_RI"], str):
                        group_ri = float(peak_group_df.loc[group, "Fiehn_RI"])
                if group_ri < RT_df["RI"].min():
                    left_idx = self.find_nearest(RT_df["RI"].values.tolist(), 0)
                    right_idx = self.find_nearest(RT_df["RI"].values.tolist(), RT_df["RI"].min())
                else:
                    left_idx = self.find_nearest(RT_df["RI"].values.tolist(), group_ri - search_wid)
                    right_idx = self.find_nearest(RT_df["RI"].values.tolist(), group_ri + search_wid)
                compound_list = RT_df.iloc[left_idx: right_idx, :].index.values.tolist()

            elif retention_score_mode == "RT":
                left_idx = self.find_nearest(RT_df["RT"].values.tolist(), group_rt_min - search_wid)
                right_idx = self.find_nearest(RT_df["RT"].values.tolist(), group_rt_min + search_wid)
                compound_list = RT_df.iloc[left_idx: right_idx, :].index.values.tolist()
            rt_index = self.find_nearest(smooth_df.index.values.tolist(), group_rt_s)
            group_rt_1 = smooth_df.index.values.tolist()[rt_index]
            sim_ion = []
            for ion in smooth_df.columns.values.tolist():
                if smooth_df.loc[group_rt_1, ion] != 0:
                    if str(ion).isdigit() == True:
                        sim_ion.append(ion)

            if peak_group_search == 1:
                compare_df = row["ions"]
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

            compare_df["intensity"] = compare_df["intensity"].apply(lambda x: (x / compare_df["intensity"].max()) * 999)
            direct_compare_df["intensity"] = direct_compare_df["intensity"].apply(
                lambda x: (x / direct_compare_df["intensity"].max()) * 999)

            if len(compound_list) > 0:
                for compound in compound_list:
                    compound_ion = []
                    if compound in meta:
                        for ion_1, intensity in meta[compound].items():
                            if isinstance(intensity, (int, float, np.float64)) or (
                                    isinstance(intensity, str) and intensity.isdigit()):
                                if int(intensity) > 0:
                                    compound_ion.append(int(ion_1))
                            else:
                                print(compound, "format error")
                    shared_ion = [i for i in compound_ion if i in sim_ion]

                    com_df_l_dic = {}
                    for i in shared_ion:
                        com_df_l_dic[i] = meta[compound][int(i)]
                    com_df_l = pd.Series(com_df_l_dic, dtype=float)
                    compare_df_1 = pd.concat([compare_df, com_df_l], axis=1, ignore_index=False)
                    compare_df_1.fillna(0, inplace=True)
                    compare_df_1 = compare_df_1.loc[~(compare_df_1 == 0).all(axis=1)]
                    compare_df_1.rename(columns={'intensity': 'intensity_s', 0: 'intensity_l'},
                                        inplace=True)
                    r_match_compare_df = compare_df_1.loc[compare_df_1.iloc[:, 1] != 0]
                    direct_compare_df_1 = pd.concat([direct_compare_df, com_df_l], axis=1, ignore_index=False)
                    direct_compare_df_1.fillna(0, inplace=True)
                    direct_compare_df_1.rename(columns={'intensity': 'intensity_s_direct', 0: 'intensity_l_direct'},
                                               inplace=True)
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
                            elif 'Fiehn_RI' in peak_group_df.columns:
                                ri_s = row['Fiehn_RI']

                            ri_l = RT_df.loc[compound, 'RI']
                            if type(ri_l) not in [float, int, np.float64, np.int64] or ri_l == np.nan:
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
                            Retention_score = 0

                    score_df.loc[compound, "Retention_score"] = Retention_score

                    if match_weight + r_match_weight == 0 or group_weight + direct_weight == 0:
                        break
                    elif MS_score * R_MS_score == 0:
                        score = 0
                    elif len(row["ions"]) == 1:
                        score = (MS_score * match_weight + R_MS_score * r_match_weight) / (
                                match_weight + r_match_weight) * 0.6 - Retention_score
                    elif len(row["ions"]) == 2:
                        score = (MS_score * match_weight + R_MS_score * r_match_weight) / (
                                match_weight + r_match_weight) * 0.7 - Retention_score
                    elif len(row["ions"]) == 3:
                        score = (MS_score * match_weight + R_MS_score * r_match_weight) / (
                                match_weight + r_match_weight) * 0.94 - Retention_score
                    elif len(row["ions"]) == 4:
                        score = (MS_score * match_weight + R_MS_score * r_match_weight) / (
                                match_weight + r_match_weight) * 0.97 - Retention_score
                    else:
                        score = (MS_score * match_weight + R_MS_score * r_match_weight) / (
                                match_weight + r_match_weight) - Retention_score

                    if score < 0:
                        score = 0

                    score_df.loc[compound, "score"] = score
                    if len(row["ions"]) >= minimum_number_of_ions and score >= sim_threshold:
                        compare_df_1 = compare_df_1.sort_index()
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
                        results_df.loc[group_rt_s, "All_match"].loc[compound, "mz_direct"] = str(
                            list(direct_compare_df_1.index))

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
            results_df.loc[group_rt_s, "All_match_list"] = [[results_list]]
            if results_df.loc[group_rt_s, "All_match"].shape[0] > 0:
                results_df.loc[group_rt_s, "All_match"].sort_values(by="Score", inplace=True, ascending=False)
                name = results_df.loc[group_rt_s, "All_match"].index.values.tolist()[0]
                results_df.loc[group_rt_s, "Best_match"].loc[name, "Score"] = \
                    results_df.loc[group_rt_s, "All_match"].iloc[
                        0, 0]
            elif results_df.loc[group_rt_s, "All_match"].shape[0] == 0 and len(row["ions"]) > 2:
                results_df.loc[group_rt_s, "Best_match_name"] = "Unknown"
            else:
                results_df.loc[group_rt_s, "Best_match_name"] = "Unknown"

            if len(score_df) > 0:
                results_df = self.New_Compound(group_rt_s, row, score_df, results_df,
                                               minimum_number_of_ions, sim_threshold)
            elif len(score_df) == 0:
                results_df.loc[group_rt_s, "New_Compound"] = ("This compound might be a new compound not present in current library")
            if len(row["ions"]) >= minimum_number_of_ions:
                results_df.loc[group_rt_s, "Relative_Peak_Area"] = peak_group_df.loc[group,
                "Relative_Peak_Area"]
                results_df.loc[group_rt_s, "Peak_Height"] = peak_group_df.loc[group, "Peak_Height"]
                results_df.loc[group_rt_s, "Quant_Ion"] = peak_group_df.loc[group, "Quant_Ion"]
        results_df_filtered = pd.DataFrame(index=results_df.index, columns=["Name", "Score"])
        for index, row1 in results_df.iterrows():
            if row1["Best_match"].shape[0] > 0:
                results_df_filtered.loc[index, "Name"] = row1["Best_match"].index.values[0]
                results_df_filtered.loc[index, "Score"] = row1["Best_match"].iloc[0, 0]
        results_df_filtered.sort_values(by='Score', inplace=True, ascending=False)
        results_df_filtered.dropna(axis=0, how="any", inplace=True)
        index_list_filter = results_df_filtered.index.values.tolist()
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

        return results_df

    def RT_to_Kovats_RI_transform(self, rt_sample, standard_df, RI_lower_limit, RI_upper_limit):
        """
        FUNCTION: Convert sample RT-RI
        rt_sample: Sample RT
        standard_df: Reference standard RT-RI
        return: Sample RI
        """
        if RI_lower_limit is None:
            RI_lower_limit = 0
        if RI_upper_limit is None:
            RI_upper_limit = 3000
        if np.isnan(rt_sample):
            return "No valid RT."
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
            z1 = np.polyfit(x, y, 5)
        return np.polyval(z1, rt_sample * 60 * 1000)

    def linear_fit(self, arg0, arg1, arg2, arg3):
        x = np.array([arg0, arg1])
        y = np.array([arg2, arg3])
        return np.polyfit(x, y, 1)

    def fill_in_missing_data(self, old_df):
        """
        Fill in the missing data errors caused by machine detection and fill them with zeros.
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
        df = pd.DataFrame(index=scan_acquisition_time, columns=list(range(30, 401)))

        if len(scan_index) == len(intensity_values) == len(mass_values) == len(scan_acquisition_time):
            r_mass_values = mass_values[::-1]
            w = 0
            while w <= len(mass_values) - 2:
                if r_mass_values[w + 1] - r_mass_values[w] < 0:
                    last_index = len(mass_values) - w - 1
                    break
                else:
                    w = w + 1

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
            rt = scan_acquisition_time[last_index]
            temp_ion = mass_values[last_index:]
            temp_intens = intensity_values[last_index:]
            for ion in temp_ion:
                ion = int(ion)
                m = np.where(temp_ion == ion)[0][0]
                df.loc[rt, ion] = temp_intens[m]

        else:
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
        df = df.dropna(axis=0, how="all")
        df = df.fillna(0)
        df = self.fill_in_missing_data(df)
        df.index = df.index.astype(float)
        return df

    def read_data_mzML(self, filename):
        run = pymzml.run.Reader(filename)
        df = pd.DataFrame()
        for spectrum in run:
            rt = (spectrum.scan_time_in_minutes()) * 60
            mz_int_dict = {}
            for mz, intensity in spectrum.peaks('centroided'):
                mz_int_dict[int(mz)] = intensity
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
        RI_min = 0
        qualitative_and_quantitative_analysis_result = pd.DataFrame()
        smooth_df = pd.DataFrame()
        total_result_file = total_result_pkl
        total_result = None

        if os.path.exists(total_result_file):
            with open(total_result_file, "rb") as f:
                total_result = pickle.load(f)
                peak_group_df = total_result[0][4]
                qualitative_and_quantitative_analysis_result = total_result[1]
                smooth_df = total_result[0][0]

        if total_result is None:
            peak_group_result = None
            if peak_group_pkl != '':
                smooth_df, smooth_df_final, decon_peak_dic, decon_data_df, peak_group_df, quant_result_df, df = peak_group_pkl
                peak_group_result = peak_group_pkl

            else:
                if isinstance(filename, pd.DataFrame):
                    df = filename
                else:
                    if filename.endswith(".mzML"):
                        df = self.read_data_mzML(filename)
                    elif filename.endswith(".cdf"):
                        data = nc.Dataset(filename)
                        df = self.read_data(data)
                RT_max = max(df.index) / 60
                smooth_df = self.lwma(df, smooth_value)
                smooth_df, noise_df = self.derivative(smooth_df, 10)
                ion_list = df.columns
                #ion_list = [110,111,238]
                peak_dic, con_peak_dic, smooth_df_final = self.find_peak(ion_list, smooth_df, noise_df,
                                                                         peak_filt_value, df)
                decon_peak_dic, decon_data_df = self.decon(ion_list, smooth_df_final, con_peak_dic)
                peak_dic, decon_peak_dic = self.calculate_sv(smooth_df_final, peak_dic, decon_data_df, decon_peak_dic)
                matched_wave, peak_group_df, quant_result_df = self.group_peak(peak_dic, decon_peak_dic, smooth_df,
                                                                               smooth_df_final,
                                                                               decon_data_df, 5,
                                                                               0.5, group_peak_factor, 0)
                del peak_dic, noise_df

                peak_group_result = smooth_df, smooth_df_final, decon_peak_dic, decon_data_df, peak_group_df, quant_result_df, df
            gc.collect()
            qualitative_and_quantitative_analysis_filename = "qualitative_and_quantitative_analysis.pkl"
            qualitative_and_quantitative_analysis_result = None

            if os.path.exists(qualitative_and_quantitative_analysis_filename):
                with open(qualitative_and_quantitative_analysis_filename, "rb") as f:
                    qualitative_and_quantitative_analysis_result = pickle.load(f)
            if qualitative_and_quantitative_analysis_result is None:
                meta, RI_df = self.read_msp(msp)
                RT_df_filename = ''
                RI_data_filename = ''
                if chooseinfo == "RT":
                    RT_df_filename = RTRIinfo
                elif chooseinfo == "RI":
                    RI_data_filename = RTRIinfo
                if RT_df_filename != '':
                    try:
                        RT_df = pd.read_csv(RT_df_filename, sep=",", index_col=0)
                    except UnicodeDecodeError:
                        try:
                            RT_df = pd.read_csv(RT_df_filename, sep=",", index_col=0, encoding='gbk')
                        except UnicodeDecodeError as e:
                            print("Error:", e)

                    duplicated_index = RT_df[RT_df.index.duplicated()].index

                    error_df = RT_df.loc[duplicated_index, :].copy()

                    RT_df = RT_df[~RT_df.index.duplicated(keep='first')]
                    RT_df = RT_df.sort_values(by='RT', ascending=True)
                    search_wid = page3RT_search_wid
                    match_weight = page3RT_match_weight
                    r_match_weight = page3RT_r_match_weight
                    group_weight = page3RT_group_weight
                    direct_weight = page3RT_direct_weight
                    minimum_number_of_ions = page3RT_group_minimum_number_of_ions
                    rt_participate = page3RT_ri_participate
                    sim_threshold = page3RT_sim_threshold
                    if rt_participate == True:
                        window = page3RT_window
                        level_factor = page3RT_level_factor
                        max_penalty = page3RT_max_penalty
                        no_info_penalty = page3RT_no_info_penalty
                    if rt_participate == False:
                        window = 0.3
                        level_factor = 0.05
                        max_penalty = 0
                        no_info_penalty = 0

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

                elif RI_data_filename != '':
                    try:
                        standard_df = pd.read_csv(RI_data_filename, sep=",")
                    except UnicodeDecodeError:
                        try:
                            standard_df = pd.read_csv(RI_data_filename, sep=",", encoding='gbk')
                        except UnicodeDecodeError as e:
                            print("Error:", e)
                    RI_max = page3RI_RI_max

                    peak_group_df = self.peak_group_df_add_retention_infor(peak_group_df, "Kovats", standard_df, RI_min,
                                                                           RI_max)

                    duplicated_index = RI_df[RI_df.index.duplicated()].index

                    error_df = RI_df.loc[duplicated_index, :].copy()

                    RI_df = RI_df[~RI_df.index.duplicated(keep='first')]
                    RI_df = RI_df.sort_values(by='RI', ascending=True)
                    search_wid = page3RI_search_wid
                    match_weight = page3RI_match_weight
                    r_match_weight = page3RI_r_match_weight
                    group_weight = page3RI_group_weight
                    direct_weight = page3RI_direct_weight
                    minimum_number_of_ions = page3RI_group_minimum_number_of_ions
                    ri_participate = page3RI_ri_participate
                    sim_threshold = page3RI_sim_threshold
                    if ri_participate == True:
                        window = page3RI_window
                        window_scale = page3RI_window_scale
                        level_factor = page3RI_level_factor
                        max_penalty = page3RI_max_penalty
                        no_info_penalty = page3RI_no_info_penalty
                        inaccurate_ri_threshold = page3RI_inaccurate_ri_threshold
                        inaccurate_ri_level_factor = page3RI_inaccurate_ri_level_factor

                    if ri_participate == False:
                        window = 100
                        window_scale = 0
                        level_factor = 0
                        max_penalty = 0
                        no_info_penalty = 0
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
                else:
                    match_weight = page3none_match_weight
                    r_match_weight = page3none_r_match_weight
                    group_weight = page3none_group_weight
                    direct_weight = page3none_direct_weight
                    minimum_number_of_ions = page3none_group_minimum_number_of_ions
                    sim_threshold = page3none_sim_threshold
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
            total_result = (peak_group_result, qualitative_and_quantitative_analysis_result)
        total_result_pkl = pickle.dumps(total_result)
        return total_result_pkl