import os
import re
import pandas as pd
import numpy as np


class CombineRtMsp():

    def combine_msp_file(self, msp_path, out_path):
        """
        Combine multiple MSP files into one combined MSP file.

        Args:
            msp_path (list): List of MSP file paths to be merged.
            out_path (str): Output directory path for the combined MSP file.
        """
        file_extension = '/*.msp'
        with open(out_path + '/combine_data.msp', 'w+') as outfile:
            for filename in msp_path:
                with open(filename) as infile:
                    outfile.write(infile.read())
        outfile.close()

    def group_cmp_inf(self, lines):
        """
        Retrieve the lines where compound names are found in an MSP file.

        Args:
            lines (list): List of lines from the MSP file.

        Returns:
            list: List of indices indicating the lines where compound names are found.

        """
        group_inf_idx = [i for i, p in enumerate(lines) if 'Name:' in p]
        group_inf_idx.append(len(lines))
        return group_inf_idx

    def del_none_ion_cpm(self, msp, error_df):
        """
        Remove compounds without ion information from an MSP file.

        Args:
            msp (file): The MSP file to process.
            error_df (pd.DataFrame): A DataFrame to store error information.

        Returns:
            list: Updated list of lines from the MSP file.
            pd.DataFrame: Updated error information DataFrame.

        """
        lines = msp.readlines()
        new_list = [item.replace('NAME:', 'Name:') for item in lines]
        lines = [item.replace("Num peaks:", "Num Peaks:") for item in new_list]
        group_inf_idx = self.group_cmp_inf(lines)
        del_list = []
        del_none_ion_cpm_list = []
        for j in range(len(group_inf_idx) - 1):
            group_inf = lines[group_inf_idx[j]:group_inf_idx[j + 1]]
            result = [item for item in group_inf if re.match(r'^\d', item)]
            if not result:
                del_none_ion_cpm_list.append(lines[group_inf_idx[j]])
                error_df.loc[len(error_df.index)] = [lines[group_inf_idx[j]].replace('Name: ', ''),
                                                     'WARNING: The compound was not found in the provided MSP file.']
                del_list.extend(
                    (
                        group_inf_idx[j],
                        group_inf_idx[j + 1],
                    )
                )

        del_com = len(del_list) - 1
        while del_com > 0:
            del lines[del_list[del_com - 1]:del_list[del_com]]
            del_com = del_com - 2
        return lines, error_df

    def Replace_Greek_numbers(self, lines):
        """
        Replace Greek letters in lines with their equivalent names.

        Args:
            lines (list): List of lines to search and replace Greek letters.

        Returns:
            list: Updated list of lines with Greek letters replaced by names.

        """
        replacements = {".alpha.": "alpha", ".beta.": "beta", ".gamma.": "gamma", ".delta.": "delta",
                        ".omega.": "omega", ".tau.": "tau"}

        for i in range(len(lines)):
            for old_str, new_str in replacements.items():
                lines[i] = lines[i].replace(old_str, new_str)
        return lines

    def find_insert_position(self, lst, num):
        left, right = 0, len(lst) - 1
        while left <= right:
            mid = (left + right) // 2
            if lst[mid] == num:
                return mid
            elif lst[mid] < num:
                left = mid + 1
            else:
                right = mid - 1
        return left

    def Remove_Duplicates(self, lines, error_df, out_path):
        """
        Remove duplicates from an MSP file.

        Args:
            lines (list): List of lines from the MSP file.
            error_df (pd.DataFrame): DataFrame to store error information.
            out_path (str): Output path for the updated MSP file.

        Returns:
            list: Updated list of lines with duplicates removed.
            pd.DataFrame: Updated error information DataFrame.

        """
        cas_strings = [s for s in lines if s.startswith('CAS')]
        cas_strings_without_dash = [s.replace('-', '') for s in cas_strings]
        cas_numbers = [int(re.findall(r'\d+', s)[0]) for s in cas_strings_without_dash]
        cas_indices = [i for i, s in enumerate(lines) if s in cas_strings]

        j_list = sorted(set(cas_numbers))
        del_list = []
        group_inf_idx = self.group_cmp_inf(lines)

        for j in j_list:
            index_list = [i for i in range(len(cas_numbers)) if cas_numbers[i] == j]
            index_list = index_list[1:]
            for idx in index_list:
                position = self.find_insert_position(group_inf_idx, cas_indices[idx])
                error_df.loc[len(error_df.index)] = [lines[group_inf_idx[position - 1]].replace('Name: ', ''),
                                                     'WARNING: This compound already exists in the database. Duplicates are not allowed.']
                del_list.extend(
                    (
                        group_inf_idx[position - 1],
                        group_inf_idx[position]
                    )
                )

        del_list.sort()
        del_com = len(del_list) - 1
        while del_com > 0:
            if del_com < len(lines):
                del lines[del_list[del_com - 1]:del_list[del_com]]
            del_com = del_com - 2
        group_inf_idx = self.group_cmp_inf(lines)
        del_list_n = []
        n_all_list = [lines[s].replace('Name: ', '') for s in group_inf_idx[:-1]]
        n_list = sorted(set(n_all_list))

        for n in n_list:
            index_list = [i for i in range(len(n_all_list)) if n_all_list[i] == n]
            index_list = index_list[1:]
            for idx in index_list:
                error_df.loc[len(error_df.index)] = [lines[group_inf_idx[idx]].replace('Name: ', ''),
                                                     'WARNING: This compound already exists in the database. Duplicates are not allowed.']
                del_list_n.extend(
                    (
                        group_inf_idx[idx],
                        group_inf_idx[idx + 1]
                    )
                )

        del_list_n.sort()
        del_com = len(del_list_n) - 1
        while del_com > 0:
            if del_com < len(lines):
                del lines[del_list_n[del_com - 1]:del_list_n[del_com]]
            del_com = del_com - 2
        with open(out_path + '/Remove_Duplicates.msp', 'w+') as f:
            for p in lines:
                f.write(p)
        f.close()
        return lines, error_df

    def combine_RT_file(self, rt_path, out_path, file_suffixes=None, merged_file_name=None):
        """
        Combine multiple RT files into a single RT library file.

        Args:
            rt_path (list): List of file paths to RT data files.
            out_path (str): Output path for the merged RT file.
            file_suffixes (list): List of file suffixes to consider for merging.
            merged_file_name (str): Name for the merged RT library file.

        Returns:
            None

        """
        merged_df = pd.DataFrame()
        if merged_file_name is None:
            merged_file_name = "combine_RT_file.xlsx"
        for file_path in rt_path:
            if file_path.split(".")[-1] == "xlsx":
                df = pd.read_excel(file_path)
                merged_df = pd.concat([merged_df, df], ignore_index=True)
            elif file_path.split(".")[-1] == "txt":
                df = pd.read_csv(file_path, sep="\t")
                merged_df = pd.concat([merged_df, df], ignore_index=True)
            elif file_path.split(".")[-1] == "csv":

                try:
                    df = pd.read_csv(file_path)
                    merged_df = pd.concat([merged_df, df], ignore_index=True)
                except UnicodeDecodeError:
                    try:
                        df = pd.read_csv(file_path, encoding='gbk')
                        merged_df = pd.concat([merged_df, df], ignore_index=True)
                    except UnicodeDecodeError as e:
                        print("Error:", e)

        merged_df.to_excel(out_path + '/' + merged_file_name, index=False)

    def Remove_Duplicates_RT(self, msp, RT_data, out_path, check_latin):
        """
        Remove duplicate entries from the RT data.

        Args:
            msp (file): MSP file containing data.
            RT_data (DataFrame): DataFrame containing RT data.
            out_path (str): Output path for error file.
            check_latin (bool): Whether to check for Latin characters.

        Returns:
            DataFrame: Cleaned RT data DataFrame.
            DataFrame: Error information DataFrame.

        """
        error_df = pd.DataFrame(columns=["Name", "RT", "Error", "NIST Preferred Name"])
        for index_1, row in RT_data.iterrows():
            try:
                float(row[1])
            except Exception:
                RT_data.drop(index=index_1, inplace=True)
                error_df.loc[len(error_df)] = [row[0], row[1], "WARNING: The RT value is out of the valid range.",
                                               np.nan]
        RT_data['RT'] = RT_data['RT'].astype('float')
        replacements = {".alpha.": "alpha", ".beta.": "beta", ".gamma.": "gamma", ".delta.": "delta",
                        ".omega.": "omega", ".tau.": "tau"}
        if check_latin:
            for i in range(RT_data.shape[0]):
                for old_str, new_str in replacements.items():
                    RT_data.iloc[i, 0] = RT_data.iloc[i, 0].replace(old_str, new_str)
        lines = msp.readlines()
        new_list = [item.replace('NAME:', 'Name:') for item in lines]
        lines = [item.replace("Num peaks:", "Num Peaks:") for item in new_list]
        lines = self.Replace_Greek_numbers(lines)
        name_df = pd.DataFrame(columns=["Name", "Index"])
        for i in lines:
            if "Name" in i:
                name_df.loc[len(name_df)] = [i.split("Name:")[1].strip(), lines.index(i)]
        name_df = name_df.drop_duplicates(subset=['Name'], keep='first')
        synon_df = pd.DataFrame(columns=["Synon", "Index"])
        for i in lines:
            if "Synon" in i:
                synon_df.loc[len(synon_df)] = [i.split("Synon:")[1].strip(), lines.index(i)]
        name_index_list = name_df["Index"].values.tolist()
        for index_1, row in RT_data.iterrows():
            if row[0] not in name_df["Name"].values.tolist():
                # print(row[0], "不在NIST首名列表")
                if row[0] not in synon_df["Synon"].values.tolist():
                    # print(row[0], "不在NIST库")
                    error_df.loc[len(error_df)] = [row[0], row[1],
                                                   "WARNING: The compound was not found in the provided MSP file.",
                                                   np.nan]
                    RT_data.drop(index=index_1, inplace=True)
                else:
                    temp_list = []
                    synon_index = synon_df.loc[synon_df["Synon"] == row[0], "Index"].tolist()[0]
                    if temp_list := [
                        x for x in name_index_list if x < synon_index
                    ]:
                        name_index = max(temp_list)
                        error_df.loc[len(error_df)] = [row[0], row[1],
                                                       "WARNING: The synonym name has been changed to NIST Preferred Name.",
                                                       name_df.loc[
                                                           name_df["Index"] == name_index, "Name"].values.tolist()[
                                                           0]]
                        RT_data.loc[index_1, "Name"] = \
                            name_df.loc[name_df["Index"] == name_index, "Name"].values.tolist()[
                                0]
                    else:
                        error_df.loc[len(error_df)] = [row[0], row[1],
                                                       "WARNING: The synonym was not found in the synonym database.",
                                                       np.nan]
                        RT_data.drop(index=index_1, inplace=True)
        duplicated_names = RT_data[RT_data.duplicated(['Name'], keep="first")]['Name']
        for name in duplicated_names:
            rt_values = RT_data.loc[RT_data['Name'] == name, 'RT'].tolist()
            name_rt_df = pd.DataFrame({'Name': [name] * len(rt_values), 'RT': rt_values,
                                       "Error": "This compound already exists in the database. Duplicates are not allowed."})
            error_df = pd.concat([error_df, name_rt_df], ignore_index=True)

        error_df.sort_values(by="Name", inplace=True, ascending=True)
        RT_data = RT_data.drop_duplicates(subset=['Name'], keep='first')
        RT_data = RT_data.sort_values(by="RT", ascending=True)
        return RT_data, error_df

    def inspection_result(self, path_rt, RT_data, standard_df, RI_alert_lower_limit, RI_alert_upper_limit,
                          RI_threshold_value, ri_window_scale, RT_lower_limit,
                          RT_upper_limit, RI_lower_limit, RI_upper_limit, check_RT):
        """
        FUNCTION: Transform RT to RI, compare RI values, and generate inspection results.

        Args:
            path_rt (str): Path to the MSP file.
            RT_data (pd.DataFrame): RT library data.
            standard_df (pd.DataFrame): Standard data.
            RI_alert_lower_limit (float): Lower limit for RI alerts.
            RI_alert_upper_limit (float): Upper limit for RI alerts.
            RI_threshold_value (float): RI threshold value.
            ri_window_scale (float): RI window scale.
            RT_lower_limit (float): Lower limit for RT.
            RT_upper_limit (float): Upper limit for RT.
            RI_lower_limit (float): Lower limit for RI.
            RI_upper_limit (float): Upper limit for RI.
            check_RT (bool): Flag to check RT.

        Returns:
            pd.DataFrame: Inspection results.
        """
        msp = open(path_rt, "r")
        lines = msp.readlines()
        new_list = [item.replace('NAME:', 'Name:') for item in lines]
        lines = [item.replace("Num peaks:", "Num Peaks:") for item in new_list]
        group_inf_idx = self.group_cmp_inf(lines)
        RI_df = pd.DataFrame(columns=['Name', 'RI_msp'])
        for j in range(len(group_inf_idx) - 1):
            group_inf = lines[group_inf_idx[j]:group_inf_idx[j + 1]]
            prefixes = [r'SemiStdNP=\d+', r'RI:\d+\n', r'Any=\d+']
            pattern = "|".join(prefixes)
            for string in group_inf:
                if 'Name:' in string:
                    RI_list = [string.replace('Name: ', '')]
                if matches := re.findall(pattern, string):
                    RI_list.extend([int(re.findall(r"\d+", match)[0]) for match in matches])
                    RI_list = RI_list[0:2]
                    RI_df.loc[len(RI_df.index)] = RI_list
        RI_df['Name'] = RI_df['Name'].str.rstrip('\n')

        combine_df = pd.merge(RT_data, RI_df, on='Name', how='outer')
        combine_df['RI_input'] = combine_df.apply(
            lambda row: self.RT_to_Kovats_RI_transform(row['RT'], standard_df, RI_lower_limit, RI_upper_limit),
            axis=1)

        combine_df['Alert'] = combine_df.apply(
            lambda row: self.alert_ri_offset_is_too_large(row['RI_msp'], row['RI_input'], RI_alert_lower_limit,
                                                          RI_alert_upper_limit, RI_threshold_value, ri_window_scale),
            axis=1)

        combine_df = combine_df[['Name', 'RT', 'RI_msp', 'RI_input', 'Alert']]

        for i, combine_df_row in combine_df.iterrows():
            if np.isnan(combine_df_row[1]):
                combine_df.loc[i, "RT"] = self.Kovats_RI_to_RT_transform(combine_df_row[2], standard_df,
                                                                         RT_lower_limit, RT_upper_limit)
                combine_df.loc[i, "Alert"] = 'rt_is_in_silico'
            if combine_df.loc[i, "Alert"] == "RI_offset_is_too_large":
                if check_RT == False:
                    combine_df.loc[i, "RT_in_silico"] = self.Kovats_RI_to_RT_transform(combine_df_row[2], standard_df,
                                                                                       RT_lower_limit,
                                                                                       RT_upper_limit)
                elif check_RT == True:
                    combine_df.loc[i, "RT"] = self.Kovats_RI_to_RT_transform(combine_df_row[2], standard_df,
                                                                             RT_lower_limit,
                                                                             RT_upper_limit)
        combine_df = combine_df.drop_duplicates(subset=['Name'], keep='first')
        combine_df = combine_df.sort_values(by="RT", ascending=True)
        combine_df.set_index(['Name'], inplace=True)
        combine_df.dropna(how='all', subset=['RT', 'RI_msp'], inplace=True)

        return combine_df

    def alert_ri_offset_is_too_large(self, RI_msp, RI_input, RI_alert_lower_limit, RI_alert_upper_limit,
                                     RI_threshold_value,
                                     ri_window_scale):
        """
        Check if the RI offset is too large.

        Args:
            RI_msp (float): RI value from the MSP file.
            RI_input (float): Input RI value.
            RI_alert_lower_limit (float): Lower limit for the RI alert.
            RI_alert_upper_limit (float): Upper limit for the RI alert.
            RI_threshold_value (float): Threshold value for the RI.
            ri_window_scale (float): RI window scale.

        Returns:
            str or None: "RI_offset_is_too_large" if the offset is too large, otherwise None.

        """
        if type(RI_input) != str:
            if RI_alert_lower_limit < RI_input < RI_alert_upper_limit and abs(
                    RI_msp - RI_input) > RI_threshold_value + ri_window_scale * 0.01 * RI_input:
                return "RI_offset_is_too_large"

    def Kovats_RI_to_RT_transform(self, ri_sample, standard_df, RT_lower_limit, RT_upper_limit):
        """
        Convert Kovats RI to RT for a given sample.

        Args:
            ri_sample (float): The RI value of the sample.
            standard_df (pd.DataFrame): Standard data containing RT-RI values.
            RT_lower_limit (float): Lower limit for RT.
            RT_upper_limit (float): Upper limit for RT.

        Returns:
            float: Transformed RT value for the sample.
        """
        prev_rows = standard_df.loc[(standard_df['RI'] < ri_sample)].tail(1)
        next_rows = standard_df.loc[(standard_df['RI'] > ri_sample)].head(1)
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
        rt_sample = RT_low + (ri_sample - RI_low) * (RT_high - RT_low) / (RI_high - RI_low)
        if rt_sample < RT_lower_limit:
            return RT_lower_limit
        elif rt_sample > RT_upper_limit:
            return RT_upper_limit
        else:
            return rt_sample

    def RT_to_Kovats_RI_transform(self, rt_sample, standard_df, RI_lower_limit, RI_upper_limit):
        """
        Convert RT to Kovats RI for a given sample.

        Args:
            rt_sample (float): The RT value of the sample.
            standard_df (pd.DataFrame): Standard data containing RT-RI values.
            RI_lower_limit (float): Lower limit for RI.
            RI_upper_limit (float): Upper limit for RI.

        Returns:
            float: Transformed Kovats RI value for the sample.
        """
        if np.isnan(rt_sample):
            return "The content of the retention time actually detected was not retrieved"
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

    def read_msp(self, msp_file):
        """
        Read and parse data from an MSP file.

        Args:
            msp_file (str): Path to the MSP file to be read.

        Returns:
            dict: A dictionary containing parsed data.
        """
        msp_file = open(msp_file, "r")
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
                    for ion in line3:
                        ion1 = ion.split(" ")
                        if len(ion1) == 2 and self.is_number(ion1[0]) and self.is_number(ion1[1]):
                            key = round(float(ion1[0]))
                            value = float(ion1[1])
                            if key in ion_intens_dic:
                                ion_intens_dic[key] = max(ion_intens_dic[key], value)
                            else:
                                ion_intens_dic[key] = value
                else:
                    print('WARNING: The file format cannot be recognized at the moment.')
            meta[name_1] = ion_intens_dic

        return meta

    def is_number(self, s):
        """
        Check if a string represents a numeric value.

        Args:
            s (str): The string to check.

        Returns:
            bool: True if the string is a valid numeric value, False otherwise.
        """
        pattern = r'^[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?$'
        return bool(re.match(pattern, s))

    def dot_product_distance(self, p, q):

        if (np.sum(p)) == 0 or (np.sum(q)) == 0:
            score = 0
        else:
            score = np.power(np.sum(q * p), 2) / \
                    (np.sum(np.power(q, 2)) * np.sum(np.power(p, 2)))
        return score

    def weighted_dot_product_distance(self, compare_df, fr_factor):
        """
        Calculate the weighted dot product distance between two spectra.

        Args:
            compare_df (pd.DataFrame): DataFrame containing spectral data.
            fr_factor (int): Factor for considering shared peaks.

        Returns:
            float: The composite score representing the distance.
        """
        m_q = pd.Series(compare_df.index)
        m_q = m_q.astype(float)
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
        if m >= fr_factor:
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

    def Main(self, msp_path, rt_path, ri_path, RI_alert_lower_limit, RI_alert_upper_limit,
             RI_threshold_value, ri_window_scale, RT_lower_limit, RT_upper_limit, RI_lower_limit, RI_upper_limit,
             check_RT, check_latin, out_path, use_unknown, unknow_msp_path, unknow_rt_path, rt_window_unknown,
             similarity_score_threshold_unknown):
        self.combine_msp_file(msp_path, out_path)

        path_msp = out_path + "/combine_data.msp"
        with open(path_msp, "r") as msp:
            error_df = pd.DataFrame(columns=['Name', 'reason'])
            lines, error_df = self.del_none_ion_cpm(msp, error_df)
            new_list = [item.replace('NAME:', 'Name:') for item in lines]
            lines = [item.replace("Num peaks:", "Num Peaks:") for item in new_list]
            if check_latin:
                lines = self.Replace_Greek_numbers(lines)
            lines, error_df = self.Remove_Duplicates(lines, error_df, out_path)
            error_df = error_df.applymap(lambda x: x.strip())
            error_df.to_excel(out_path + '/error_df.xlsx', index=False)

        self.combine_RT_file(rt_path, out_path)
        path_rt = out_path + "/Remove_Duplicates.msp"

        msp = open(path_rt, "r")
        RT_data = pd.DataFrame(columns=['Name', 'RT'])
        try:
            RT_data_file = out_path + "/combine_RT_file.xlsx"
            if os.path.exists(RT_data_file):
                RT_data = pd.read_excel(RT_data_file, header=0)
            if RT_data.empty:
                data = np.empty((0, 2))
                RT_data = pd.DataFrame(data, columns=['Name', 'RT'])
            RT_data, error_df = self.Remove_Duplicates_RT(msp, RT_data, out_path, check_latin)
            msp.close()
        except:
            msp.close()
        msp = open(path_rt, "r")

        try:
            standard_df = pd.read_csv(ri_path, sep=",")
            combine_df = self.inspection_result(path_rt, RT_data, standard_df, RI_alert_lower_limit,
                                                RI_alert_upper_limit,
                                                RI_threshold_value, ri_window_scale, RT_lower_limit, RT_upper_limit,
                                                RI_lower_limit, RI_upper_limit, check_RT)
            combine_df.to_csv(out_path + '/New_RT_list.csv', index=True)
        except:
            if rt_path != [] and ri_path != '':
                RT_data.to_csv(out_path + '/New_RT_list.csv', index=True)
            elif rt_path != [] and ri_path == '':
                RT_data.to_csv(out_path + '/New_RT_list.csv', index=True)

        if use_unknown:
            meta = self.read_msp(out_path + "/Remove_Duplicates.msp")
            unknow_rt = pd.read_csv(unknow_rt_path)
            unknow_meta = self.read_msp(unknow_msp_path)
            remove_duplicates_lines = open(out_path + "/Remove_Duplicates.msp").readlines()
            unknow_lines = open(unknow_msp_path).readlines()
            unknow_lines = [item.replace('NAME:', 'Name:') for item in unknow_lines]
            unknow_lines = [item.replace("Num peaks:", "Num Peaks:") for item in unknow_lines]
            unknow_group_inf_idx = self.group_cmp_inf(unknow_lines)
            add_unknown_list = []
            for index, row in unknow_rt.iterrows():
                name = row[0]
                rt = row[1]

                window_minute_name = combine_df[
                    (combine_df.iloc[:, 0] >= rt - rt_window_unknown) & (
                                combine_df.iloc[:, 0] <= rt + rt_window_unknown)].index
                unknow_ms_df = pd.DataFrame(list(unknow_meta[name].items()), columns=['ion', 'intensity'])
                name_index = unknow_lines.index("Name: " + name + "\n")
                next_name_index = unknow_group_inf_idx[unknow_group_inf_idx.index(name_index) + 1]
                name_lines = unknow_lines[name_index:next_name_index]
                for m in window_minute_name:
                    m_df = pd.DataFrame(list(meta[m].items()), columns=['ion', 'intensity'])
                    merged_df = pd.merge(unknow_ms_df, m_df, on='ion', how='inner')
                    score = self.weighted_dot_product_distance(merged_df, 2)
                    if score > similarity_score_threshold_unknown:
                        break
                else:
                    add_unknown_list.append(name)
                    remove_duplicates_lines = remove_duplicates_lines + name_lines

            with open(out_path + '/Remove_Duplicates.msp', 'w+') as f:
                for p in remove_duplicates_lines:
                    f.write(p)
            if add_unknown_list != []:

                path_rt = out_path + "/Remove_Duplicates.msp"
                msp = open(path_rt, "r")
                try:
                    RT_data_file = out_path + "/combine_RT_file.xlsx"
                    if os.path.exists(RT_data_file):
                        RT_data = pd.read_excel(RT_data_file, header=0)
                    if RT_data.empty:
                        data = np.empty((0, 2))
                        RT_data = pd.DataFrame(data, columns=['Name', 'RT'])
                    RT_data, error_df = self.Remove_Duplicates_RT(msp, RT_data, out_path, check_latin)

                    msp.close()
                except:
                    msp.close()
                msp = open(path_rt, "r")
                for i in add_unknown_list:
                    RT_data = RT_data.append(unknow_rt[unknow_rt["Name"] == i], ignore_index=True)
                try:
                    standard_df = pd.read_csv(ri_path, sep=",")
                    combine_df = self.inspection_result(path_rt, RT_data, standard_df, RI_alert_lower_limit,
                                                        RI_alert_upper_limit,
                                                        RI_threshold_value, ri_window_scale, RT_lower_limit,
                                                        RT_upper_limit,
                                                        RI_lower_limit, RI_upper_limit, check_RT)
                    combine_df.to_csv(out_path + '/New_RT_list.csv', index=True)
                    msp.close()
                except:
                    msp.close()
                    if rt_path != [] and ri_path != '':
                        RT_data.to_csv(out_path + '/New_RT_list.csv', index=True)
                    elif rt_path != [] and ri_path == '':
                        RT_data.to_csv(out_path + '/New_RT_list.csv', index=True)
        os.remove(out_path + "/combine_data.msp")
        os.remove(out_path + "/combine_RT_file.xlsx")
