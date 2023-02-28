result = []

def equal(anarci, dssp, i, j, result):
    if i == len(anarci) or j == len(dssp):
        return True

    if dssp[j][1] == anarci[i][1]:
        return equal(anarci, dssp, i + 1, j + 1, result)

    if dssp[j][1] == 'X':
        if equal(anarci, dssp, i + 1, j + 1, result):
            return True
        else:
            try_stay = equal(anarci, dssp, i, j + 1, result)
            if try_stay:
                result.append(j)
                return True

    return False

anarci = [('H1', 'Q'), ('H2', 'V'), ('H3', 'Q'), ('H4', 'L'), ('H5', 'V'), ('H6', 'Q'), ('H7', 'S'), ('H8', 'G'), ('H9', 'A'), ('H11', 'E'), ('H12', 'V'), ('H13', 'K'), ('H14', 'K'), ('H15', 'P'), ('H16', 'G'), ('H17', 'S'), ('H18', 'S'), ('H19', 'V'), ('H20', 'K'), ('H21', 'V'), ('H22', 'S'), ('H23', 'C'), ('H24', 'K'), ('H25', 'A'), ('H26', 'S'), ('H27', 'G'), ('H28', 'G'), ('H29', 'T'), ('H30', 'F'), ('H35', 'S'), ('H36', 'S'), ('H37', 'Y'), ('H38', 'A'), ('H39', 'I'), ('H40', 'S'), ('H41', 'W'), ('H42', 'V'), ('H43', 'R'), ('H44', 'Q'), ('H45', 'A'), ('H46', 'P'), ('H47', 'G'), ('H48', 'Q'), ('H49', 'G'), ('H50', 'L'), ('H51', 'E'), ('H52', 'W'), ('H53', 'M'), ('H54', 'G'), ('H55', 'S'), ('H56', 'I'), ('H57', 'I'), ('H58', 'P'), ('H59', 'W'), ('H62', 'F'), ('H63', 'G'), ('H64', 'T'), ('H65', 'T'), ('H66', 'N'), ('H67', 'Y'), ('H68', 'A'), ('H69', 'Q'), ('H70', 'K'), ('H71', 'F'), ('H72', 'Q'), ('H74', 'G'), ('H75', 'R'), ('H76', 'V'), ('H77', 'T'), ('H78', 'I'), ('H79', 'T'), ('H80', 'A'), ('H81', 'D'), ('H82', 'E'), ('H83', 'S'), ('H84', 'T'), ('H85', 'S'), ('H86', 'T'), ('H87', 'A'), ('H88', 'Y'), ('H89', 'M'), ('H90', 'E'), ('H91', 'L'), ('H92', 'S'), ('H93', 'S'), ('H94', 'L'), ('H95', 'R'), ('H96', 'S'), ('H97', 'E'), ('H98', 'D'), ('H99', 'T'), ('H100', 'A'), ('H101', 'V'), ('H102', 'Y'), ('H103', 'Y'), ('H104', 'C'), ('H105', 'A'), ('H106', 'R'), ('H107', 'D'), ('H108', 'S'), ('H109', 'E'), ('H113', 'Y'), ('H114', 'Y'), ('H115', 'F'), ('H116', 'D'), ('H117', 'H'), ('H118', 'W'), ('H119', 'G'), ('H120', 'Q'), ('H121', 'G'), ('H122', 'T'), ('H123', 'L'), ('H124', 'V'), ('H125', 'T'), ('H126', 'V'), ('H127', 'S'), ('H128', 'S')]

dssp = [('H1', 'X'), ('H1', 'X'), ('H2', 'V'), ('H3', 'Q'), ('H1', 'X'), ('H1', 'X'), ('H4', 'L'), ('H5', 'V'), ('H6', 'Q'), ('H7', 'S'), ('H8', 'G'), ('H9', 'X'), ('H10', 'E'), ('H11', 'V'), ('H12', 'K'), ('H13', 'K'), ('H14', 'P'), ('H15', 'G'), ('H16', 'S'), ('H17', 'S'), ('H18', 'V'), ('H19', 'K'), ('H20', 'V'), ('H21', 'S'), ('H22', 'C'), ('H23', 'K'), ('H24', 'A'), ('H25', 'S'), ('H26', 'G'), ('H27', 'G'), ('H28', 'T'), ('H29', 'F'), ('H30', 'S'), ('H31', 'S'), ('H32', 'Y'), ('H33', 'A'), ('H34', 'I'), ('H35', 'S'), ('H36', 'W'), ('H37', 'V'), ('H38', 'R'), ('H39', 'Q'), ('H40', 'A'), ('H41', 'P'), ('H42', 'G'), ('H43', 'Q'), ('H44', 'G'), ('H45', 'L'), ('H46', 'E'), ('H47', 'W'), ('H48', 'M'), ('H49', 'G'), ('H50', 'S'), ('H51', 'I'), ('H52', 'I'), ('H53', 'P'), ('H54', 'W'), ('H55', 'F'), ('H56', 'G'), ('H57', 'T'), ('H58', 'T'), ('H59', 'N'), ('H60', 'Y'), ('H61', 'A'), ('H62', 'Q'), ('H63', 'K'), ('H64', 'F'), ('H65', 'Q'), ('H66', 'G'), ('H67', 'R'), ('H68', 'V'), ('H69', 'T'), ('H70', 'I'), ('H71', 'T'), ('H72', 'A'), ('H73', 'D'), ('H74', 'E'), ('H75', 'S'), ('H76', 'T'), ('H77', 'S'), ('H78', 'T'), ('H79', 'A'), ('H80', 'Y'), ('H81', 'M'), ('H82', 'E'), ('H83', 'L'), ('H84', 'S'), ('H85', 'S'), ('H86', 'L'), ('H87', 'R'), ('H88', 'S'), ('H89', 'E'), ('H90', 'D'), ('H91', 'T'), ('H92', 'A'), ('H93', 'V'), ('H94', 'Y'), ('H95', 'Y'), ('H96', 'C'), ('H97', 'A'), ('H98', 'R'), ('H99', 'D'), ('H100', 'S'), ('H101', 'E'), ('H102', 'Y'), ('H103', 'Y'), ('H104', 'F'), ('H105', 'D'), ('H106', 'H'), ('H107', 'W'), ('H108', 'G'), ('H109', 'Q'), ('H110', 'G'), ('H111', 'T'), ('H112', 'L'), ('H113', 'V'), ('H114', 'T'), ('H115', 'V'), ('H116', 'S'), ('H117', 'S'), ('H118', 'A'), ('H119', 'S'), ('H120', 'T'), ('H121', 'K'), ('H122', 'G'), ('H123', 'P'), ('H124', 'S'), ('H125', 'V'), ('H126', 'F'), ('H127', 'P'), ('H128', 'L'), ('H129', 'A'), ('H130', 'P'), ('H131', 'S'), ('H132', 'S'), ('H133', 'K'), ('H134', 'S'), ('H135', 'T'), ('H136', 'S'), ('H137', 'G'), ('H138', 'G'), ('H139', 'T'), ('H140', 'A'), ('H141', 'A'), ('H142', 'L'), ('H143', 'G'), ('H144', 'C'), ('H145', 'L'), ('H146', 'V'), ('H147', 'K'), ('H148', 'D'), ('H149', 'Y'), ('H150', 'F'), ('H151', 'P'), ('H152', 'E'), ('H153', 'P'), ('H154', 'V'), ('H155', 'T'), ('H156', 'V'), ('H157', 'S'), ('H158', 'W'), ('H159', 'N'), ('H160', 'S'), ('H161', 'G'), ('H162', 'A'), ('H163', 'L'), ('H164', 'T'), ('H165', 'S'), ('H166', 'G'), ('H167', 'V'), ('H168', 'H'), ('H169', 'T'), ('H170', 'F'), ('H171', 'P'), ('H172', 'A'), ('H173', 'V'), ('H174', 'L'), ('H175', 'Q'), ('H176', 'S'), ('H177', 'S'), ('H178', 'G'), ('H179', 'L'), ('H180', 'Y'), ('H181', 'S'), ('H182', 'L'), ('H183', 'S'), ('H184', 'S'), ('H185', 'V'), ('H186', 'V'), ('H187', 'T'), ('H188', 'V'), ('H189', 'P'), ('H190', 'S'), ('H191', 'S'), ('H192', 'S'), ('H193', 'L'), ('H194', 'G'), ('H195', 'T'), ('H196', 'Q'), ('H197', 'T'), ('H198', 'Y'), ('H199', 'I'), ('H200', 'C'), ('H201', 'N'), ('H202', 'V'), ('H203', 'N'), ('H204', 'H'), ('H205', 'K'), ('H206', 'P'), ('H207', 'S'), ('H208', 'N'), ('H209', 'T'), ('H210', 'K'), ('H211', 'V'), ('H212', 'D'), ('H213', 'K'), ('H214', 'K'), ('H215', 'V'), ('H216', 'E'), ('H217', 'P')]

print(equal(anarci, dssp, 0, 0, result))

print(result)