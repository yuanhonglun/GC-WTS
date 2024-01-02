import pandas as pd
import numpy as np
import re
import os


class GetMethod():

    def replace_name(self, list_1):
        """
        Replace specific strings in a list with their corresponding replacements.

        Args:
            list_1 (list): A list of strings to be processed.

        Returns:
            list: A list with the specified replacements applied.

        """
        replacements = {".alpha.": "alpha", ".beta.": "beta", ".gamma.": "gamma", ".delta.": "delta",
                        ".omega.": "omega", ".tau.": "tau"}
        for i in range(len(list_1)):
            for old_str, new_str in replacements.items():
                list_1[i] = list_1[i].replace(old_str, new_str)

        return list_1

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

    def remove_content_from_list(self, lst, content):
        """
        Remove specific content from a list of strings using regular expressions.

        Args:
            lst (list): A list of strings to be processed.
            content (str): The content to be removed from each string in the list.

        Returns:
            list: A list with the specified content removed from the strings.

        """
        pattern = re.escape(content)
        regex = re.compile(pattern)
        return [regex.sub('', string) for string in lst if isinstance(string, str)]

    def read_msp(self, msp_file):
        """
        Read data from an MSP file and convert it into a dictionary format.

        Args:
            msp_file (str): The path to the MSP file.

        Returns:
            dict: A dictionary where keys are compound names and values are dictionaries of ion intensities.

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

        return meta

    def dot_product_distance(self, p, q):
        """
        Calculate the dot product distance between two vectors p and q.

        Args:
            p (numpy.ndarray): First vector.
            q (numpy.ndarray): Second vector.

        Returns:
            float: Dot product distance between the two vectors.
        """

        if (np.sum(p)) == 0 or (np.sum(q)) == 0:
            score = 0
        else:
            score = np.power(np.sum(q * p), 2) / \
                    (np.sum(np.power(q, 2)) * np.sum(np.power(p, 2)))
        return score

    def weighted_dot_product_distance(self, compare_df, fr_factor):

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

    def calculate_similarity(self, target_name, df, n, fr_factor):
        """
        Calculate the similarity scores between the target compound and all compounds in the DataFrame.

        Args:
            target_name (str): The name of the target compound.
            df (pd.DataFrame): The DataFrame containing compound data.
            n (int): The number of top similar compounds to return.
            fr_factor (float): The factor used in the weighted dot product distance calculation.

        Returns:
            pd.DataFrame: A DataFrame containing similarity scores for each compound.

        """
        result_df = pd.DataFrame(columns=["Score"])
        first_col = df.loc[target_name]
        for compound in df.index.values:
            if compound != target_name:
                second_col = df.loc[compound]
                compare_df = pd.concat([first_col, second_col], axis=1)
                compare_df = compare_df.astype(float)
                score = self.weighted_dot_product_distance(compare_df, fr_factor)
                result_df.loc[compound, "Score"] = score

        return result_df

    def calculate_average_score_and_difference_count(self, targeted_compound,
                                                     ion_combination,
                                                     df, similarity_threshold,
                                                     fr_factor,
                                                     n):
        """
        Calculate the average similarity score and the difference count for a targeted compound using specified ion combinations.

        Args:
            targeted_compound (str): The name of the targeted compound.
            ion_combination (list): List of ion combinations for similarity calculation.
            df (pd.DataFrame): The DataFrame containing compound data.
            similarity_threshold (float): The similarity threshold for considering compounds as similar.
            fr_factor (float): The factor used in the weighted dot product distance calculation.
            n (int): The number of top similar compounds to consider.

        Returns:
            pd.DataFrame: A DataFrame containing difference counts and average similarity scores for each ion combination.

        """
        difference_count_df = pd.DataFrame(columns=["Diff_Count", "Similar_Compound_Ave_Score"])
        for ions in ion_combination:
            temp_df_1 = df[ions]
            result_df_2 = self.calculate_similarity(targeted_compound, temp_df_1, n, fr_factor)
            result_df_3 = result_df_2[(result_df_2["Score"] < similarity_threshold)]
            count = len(result_df_3)
            difference_count_df.loc[str(ions), "Diff_Count"] = count

            result_df_4 = result_df_2[(result_df_2["Score"] >= similarity_threshold)]
            if result_df_4.shape[0] > 0:
                ave_score_1 = np.average(result_df_4, axis=0)[0]
            else:
                ave_score_1 = 1
            difference_count_df.loc[str(ions), "Similar_Compound_Ave_Score"] = ave_score_1

        difference_count_df.sort_values(by="Diff_Count", inplace=True, ascending=False)
        return difference_count_df

    def calculate_combination_score(self, combination_df, targeted_compound, temp_df, prefer_mz_threshold):
        """
        Calculate the combination score for ion combinations in a DataFrame.

        Args:
            combination_df (pd.DataFrame): DataFrame containing ion combinations.
            targeted_compound (str): The name of the targeted compound.
            temp_df (pd.DataFrame): Temporary DataFrame containing compound data.
            prefer_mz_threshold (int): The preferred m/z threshold.

        Returns:
            pd.DataFrame: DataFrame with combination scores added.

        """
        for index, row in combination_df.iterrows():
            ion_list = re.findall("\d+\.?\d*", index)
            ion_list = list(map(int, ion_list))
            new_temp_df = temp_df.loc[str(targeted_compound), ion_list].to_frame()
            new_temp_df["ion"] = new_temp_df.index.tolist()
            new_temp_df["ion"] = new_temp_df["ion"].astype('int')
            new_temp_df["ion"] = np.where(new_temp_df["ion"] < prefer_mz_threshold, 1, new_temp_df["ion"])
            new_temp_df["score"] = (pow(new_temp_df["ion"], 3)) * (pow(new_temp_df[str(targeted_compound)], 0.5))
            combination_df.loc[index, "com_score"] = new_temp_df["score"].sum()

        return combination_df

    def calculate_solo_compound_combination_score(self, matrix_1, prefer_mz_threshold):
        """
        Calculate the combination score for solo compounds in a DataFrame.

        Args:
            matrix_1 (pd.DataFrame): DataFrame containing solo compound data.
            prefer_mz_threshold (int): The preferred m/z threshold.

        Returns:
            pd.DataFrame: DataFrame with combination scores added and sorted by score.

        """
        matrix_1['ion'] = matrix_1['ion'].apply(lambda x: 1 if x < prefer_mz_threshold else x)
        matrix_1['com_score'] = matrix_1.apply(lambda row: pow(row.iloc[0], 0.5) * pow(row.iloc[1], 3), axis=1)
        matrix_1 = matrix_1.sort_values(by='com_score', ascending=False)
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

    def Main(self, msp_path, rt_data_path, set_name_list, name_list_path, mz_min, mz_max, outpath, rt_window,
             min_ion_intensity_percent, min_ion_num, prefer_mz_threshold, similarity_threshold, fr_factor,
             retention_time_max, solvent_delay, sim_sig_max, min_dwell_time,
             point_per_s, convert_to_ag_method):

        msp = msp_path

        try:
            RT_data = pd.read_csv(rt_data_path, index_col=0)
        except UnicodeDecodeError:
            try:
                RT_data = pd.read_csv(rt_data_path, index_col=0, encoding='gbk')
            except UnicodeDecodeError as e:
                print("Error:", e)

        error_df = pd.DataFrame(columns=["error"])
        meta_1 = self.read_msp(msp)

        matrix = pd.DataFrame()
        for name, dic in meta_1.items():
            b = {name: dic}
            y = pd.DataFrame.from_dict(b).T
            matrix = pd.concat([matrix, y], axis=0, join='outer').fillna(0)

        for col_name in list(matrix):
            if type(col_name) not in [float, int, np.float64, np.int64] and col_name.isdigit() == False:
                matrix.drop(columns=col_name, inplace=True)

        for ion in list(matrix):
            if int(ion) < mz_min or int(ion) > mz_max:
                matrix.drop(columns=ion, inplace=(True))

        matrix_name = matrix.index.tolist()

        duplicated_index = RT_data.index[RT_data.index.duplicated()]
        for index in duplicated_index:
            error_df.loc[index, "error"] = "Duplicated RT"
        RT_data.drop(duplicated_index, axis=0, inplace=True)
        for index_1, row in RT_data.iterrows():

            if index_1 not in matrix_name:
                error_df.loc[index_1, "error"] = "This compound is in RT list, but not in MSP library"
                RT_data.drop(index=index_1, inplace=True)

            elif type(row[0]) not in [float, int, np.float64, np.int64] or row[0] == np.nan:
                error_df.loc[index_1, "error"] = "RT format error"
                RT_data.drop(index=index_1, inplace=True)

        RT_data = RT_data.sort_values(by='RT')

        if set_name_list:
            f1 = open(name_list_path)
            f2 = list(f1)
            compound_list = []
            for name in f2:
                name = name.strip()
                name = name.strip('"')
                compound_list.append(name)
            compound_list = list(set(compound_list))
            for i in compound_list:
                if i not in RT_data.index.values:
                    error_df.loc[i, "error"] = "This compound is in name list, but not in RT list"
                    compound_list.remove(i)

        else:
            compound_list = RT_data.index.values.tolist()
        error_df.to_csv(outpath + "/{}.csv".format("input_data_error_info"), index=True, index_label="Name")
        nearby_compound_dic = {}
        for name in compound_list:
            if name in RT_data.index.values.tolist():
                rt = RT_data.at[name, "RT"]
                nearby_compound_dic[name] = RT_data[(RT_data.iloc[:, 0] >= rt - rt_window) &
                                                    (RT_data.iloc[:, 0] <= rt + rt_window)].index.tolist()
        combination_result_df_file = "combination_results.csv"
        combination_result_df = None

        if os.path.exists(combination_result_df_file):
            combination_result_df = pd.read_csv(r"./combination_results.csv", header=0, index_col=0)

        if combination_result_df is None:
            combination_result_df = pd.DataFrame(
                columns=["RT", "Ion_Combination", "Note", "Similar_Compound_List", "SCL_Note"])

            min_ion_intensity = min_ion_intensity_percent * 10

            for targeted_compound, nearby_compound_list in nearby_compound_dic.items():
                combination_result_df.loc[targeted_compound, "RT"] = RT_data.loc[targeted_compound, "RT"]
                title = self.validateTitle(str(targeted_compound))
                solo_list = []
                error_list_2 = []
                if nearby_compound_list == [targeted_compound]:
                    scl = []
                    combination_result_df.loc[targeted_compound, "Similar_Compound_List"] = scl
                    combination_result_df.loc[targeted_compound, "SCL_Note"] = "No adjacent compounds."

                    matrix_1 = matrix.loc[targeted_compound].to_frame()
                    matrix_1["ion"] = matrix_1.index.tolist()
                    matrix_1["ion"] = matrix_1["ion"].astype(int)
                    matrix_1[targeted_compound] = matrix_1[targeted_compound].astype(float)
                    matrix_1[targeted_compound] = np.where(matrix_1[targeted_compound] < min_ion_intensity, 0,
                                                           matrix_1[targeted_compound])
                    matrix_1 = matrix_1.loc[matrix_1[targeted_compound] > 0, :]

                    if matrix_1.shape[0] < 2:
                        combination_result_df.loc[targeted_compound, "Ion_Combination"] = "NA"
                        combination_result_df.loc[targeted_compound, "Note"] = ("No adjacent compounds; "
                                                                                "The available number of ions is less "
                                                                                "than 2, the compound is excluded")
                    else:
                        matrix_1 = self.calculate_solo_compound_combination_score(matrix_1, prefer_mz_threshold)

                        if matrix_1.shape[0] <= min_ion_num:
                            combination_list = matrix_1.index.values.tolist()
                        else:
                            combination_list = matrix_1.iloc[0: min_ion_num, :].index.values.tolist()
                        combination_result_df.loc[str(targeted_compound), "Ion_Combination"] = combination_list
                        solo_list.append(targeted_compound)
                else:
                    temp_df = matrix.loc[nearby_compound_list]
                    temp_df = temp_df.astype(float)
                    temp_df.loc[targeted_compound, :] = np.where(temp_df.loc[targeted_compound, :] < min_ion_intensity,
                                                                 0,
                                                                 temp_df.loc[targeted_compound, :])
                    temp_df = temp_df.loc[:, temp_df.loc[targeted_compound, :] > 0]

                    if temp_df.shape[1] < 2:
                        combination_result_df.loc[targeted_compound, "Ion_Combination"] = "NA"
                        combination_result_df.loc[targeted_compound, "Note"] = ("The available number of ions is less "
                                                                                "than 2, the compound is excluded")
                    else:
                        similar_compound_list = []
                        result_df_1 = self.calculate_similarity(targeted_compound, temp_df, -1, fr_factor)
                        for index, row in result_df_1.iterrows():
                            if float(row) >= similarity_threshold:
                                similar_compound_list.append(index)

                        combination_result_df.loc[targeted_compound, "Similar_Compound_List"] = similar_compound_list
                        temp_df.drop(index=similar_compound_list, inplace=True)

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
                                    combination_result_df.loc[targeted_compound, "Ion_Combination"] = "NA"
                                    combination_result_df.loc[targeted_compound, "Note"] = ("The available number of ions is less "
                                                                                "than 2, the compound is discarded")
                                else:
                                    matrix_1 = self.calculate_solo_compound_combination_score(matrix_1,
                                                                                              prefer_mz_threshold)

                                    if matrix_1.shape[0] <= min_ion_num:
                                        combination_list = matrix_1.index.values.tolist()
                                    else:
                                        combination_list = matrix_1.iloc[0: min_ion_num, :].index.values.tolist()
                                    combination_result_df.loc[
                                        str(targeted_compound), "Ion_Combination"] = combination_list
                                    solo_list.append(targeted_compound)
                        else:
                            col_name = list(temp_df.columns)
                            col_name = list(map(int, col_name))
                            new_com = [[x] for x in col_name]

                            difference_count_df_1 = self.calculate_average_score_and_difference_count(targeted_compound,
                                                                                                      new_com,
                                                                                                      temp_df,
                                                                                                      similarity_threshold,
                                                                                                      fr_factor,
                                                                                                      0)

                            combination_df = difference_count_df_1[difference_count_df_1["Diff_Count"] >=
                                                                   difference_count_df_1.iat[0, 0]]
                            if combination_df.shape[0] > 5:
                                combination_df = combination_df.sort_values(by="Similar_Compound_Ave_Score",
                                                                            inplace=False, ascending=True)
                                combination_df = combination_df[combination_df["Similar_Compound_Ave_Score"] >=
                                                                combination_df.iat[0, 1]]
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
                            ion_list = list(temp_df)
                            combination_array = combination_df.index.values
                            n = 0
                            ion_num = 1

                            flag = True

                            while True:
                                if int((combination_df.max())[0]) >= int(
                                        temp_df.shape[0] - 1) and ion_num >= min_ion_num:
                                    break
                                elif flag == False:
                                    break
                                else:
                                    n = n + 1
                                    total_list = []
                                    new_total = []

                                    for ion_combination in combination_array:

                                        ion_combination_list = re.findall("\d+\.?\d*", ion_combination)
                                        ion_combination_list = list(map(int, ion_combination_list))
                                        candidate_list = [i for i in ion_list if i not in ion_combination_list]

                                        if candidate_list == []:
                                            if int((combination_df.max())[0]) >= int(temp_df.shape[0] - 1):
                                                if combination_df.shape[0] > 1:
                                                    combination_df = combination_df.sort_values(
                                                        by="Similar_Compound_Ave_Score", inplace=False, ascending=True)
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
                                                    str(targeted_compound), "Note"] = ("Despite the qualitative ion number is less than the defined number, its separation score reaches 1")
                                                flag = False

                                                break

                                            else:
                                                combination_result_df.loc[
                                                    str(targeted_compound), "Ion_Combination"] = "NA"
                                                combination_result_df.loc[
                                                    str(targeted_compound), "Note"] = ("No ions available, this compound is discarded")
                                                flag = False
                                                break

                                        elif flag == True:

                                            for candidate in candidate_list:
                                                temp_ion_combination_list = ion_combination_list.copy()
                                                temp_ion_combination_list.append(candidate)
                                                total_list.append(temp_ion_combination_list)
                                                temp_ion_combination_list = []
                                            for i in total_list:
                                                i = list(map(int, i))
                                                i.sort()
                                                if i not in new_total:
                                                    new_total.append(i)
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
                                            combination_result_df.loc[str(targeted_compound), "Ion_Combination"] = "NA"
                                            combination_result_df.loc[
                                                str(targeted_compound), "Note"] = ("Error: The 'difference_count_df' is "
                                                                                   "empty.")
                                            flag = False

                                            break

                                        if combination_df.shape[0] > 1:
                                            combination_df = combination_df.sort_values(by="Similar_Compound_Ave_Score",
                                                                                        inplace=False, ascending=True)
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
                                        combination_array = combination_df.index.values
                                        ion_num = ion_num + 1

                                    else:

                                        break
                            if flag == True:
                                combination_result_df.loc[str(targeted_compound), "Ion_Combination"] = \
                                    combination_array[0]
            combination_result_df.to_csv(outpath + "/combination_results.csv", index=True, index_label="Name")

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
                    error_df.loc[len(error_df)] = [name, "The ion group format is incorrect."]
            else:
                error_df.loc[len(error_df)] = [name, "This compound is not in the RT list."]

        length = len(name_list_total)

        data = {'RT': RT_list_total, 'ion': num}

        ion_rt = pd.DataFrame(data, index=name_list_total)

        ion_rt.sort_values(by="RT", inplace=True, ascending=True)
        for idx, row in ion_rt.iterrows():
            if row["RT"] > retention_time_max:
                ion_rt.loc[idx, "RT"] = retention_time_max

        ion_rt.to_csv(outpath + "/ion_rt_data.csv", index=True, index_label="Name")
        rt_index = [i * 0.5 / 60 for i in range(0, int(retention_time_max * 120) + 1, 1)]
        df = pd.DataFrame(index=rt_index, columns=[i for i in range(mz_min, mz_max + 1)])
        for i, row in ion_rt.iterrows():
            df.loc[(df.index > row[0] - rt_window) & (df.index < row[0] + rt_window), row[1]] = 1

        df = df[df.index > solvent_delay]

        df["sum"] = df.sum(axis=1)
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
                        print(
                            "An error occurred when attempting to merge with 0 data points, and it resulted in the "
                            "following error message: group_id =",
                            i)

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
            if result_df.shape[0] > 1:
                result_df.sort_values(by="row_number", inplace=True, ascending=False)
                result_df = result_df[:1]
            else:
                result_df = result_df
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

        df = df[df['sum'] != 0]

        df["dwell_time"] = ""

        for i in df['group_id'].values.tolist():

            ion_num = df.loc[df["group_id"] == i, "sum"].values[0]
            if (1000 / point_per_s) / ion_num >= min_dwell_time:

                df.loc[df["group_id"] == i, "dwell_time"] = (1000 / point_per_s) / ion_num

            else:

                df.loc[df["group_id"] == i, "dwell_time"] = min_dwell_time

        df.to_csv(outpath + '/SIM_seg_result.csv', index=True)
        if convert_to_ag_method:

            from lxml import etree

            group_list = df['group_id'].unique().tolist()
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

            with open(outpath + "/qqqacqmethod.xml", "wb") as f:
                f.write(etree.tostring(root, encoding="UTF-8", xml_declaration=True))
