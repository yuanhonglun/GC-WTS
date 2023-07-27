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

def test_push():
    pass

def replace_1(self, x):
    if 1 in x.values:
        x[:] = 1
    return x
