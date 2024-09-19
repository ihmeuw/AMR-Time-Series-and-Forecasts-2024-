from builtins import range

import pandas as pd

from cod_prep.downloaders.ages import get_ages
from cod_prep.utils import print_log_message, report_if_merge_fail
from mcod_prep.datasets.mcod_formatting import MultipleCauseFormatter, validate_cause_column

DATA_COL = ["deaths"]

VALUE_COLS = ["deaths"]

SOURCE = "SOURCE"

YEARS = list(range(1998, 2023))


def get_file_path(year):
    """Get file path based on year."""
    if year in range(1998, 2012):
        date = "DATE"
    elif year in range(2012, 2016):
        date = "DATE"
    elif year in range(2016, 2018):
        date = "DATE"
    elif year == 2018:
        date = "DATE"
    elif year == 2019:
        date = "DATE"
    elif year == 2020:
        date = "DATE"
    elif year == 2021:
        date = "DATE"
    elif year == 2022:
        date = "DATE"

    if year < 2018:
        file_path = "FILEPATH"
    elif year == 2018:
        file_path = "FILEPATH"
    elif year in [2019, 2020, 2021]:
        file_path = "FILEPATH"
    else:
        file_path = "FILEPATH"

    return file_path


def get_columns(year):
    if year < 2019:
        columns = [
            "gru_ed1",
            "sexo",
            "sit_defun",
            "c_dir1",
            "c_ant1",
            "c_ant2",
            "c_ant3",
            "c_pat1",
            "c_bas1",
        ]

        if year in range(2008, 2019):
            columns += ["c_dir12", "c_ant12", "c_ant22", "c_ant32", "c_pat2"]

        if year in list(range(2002, 2004)) + list(range(2008, 2014)) + list(range(2015, 2018)):
            columns = [col.upper() if col != "sexo" else col for col in columns]

        if year == 2018:
            columns = [col.upper() for col in columns]
    else:
        columns = ["GRU_ED1", "SEXO", "SIT_DEFUN", "C_BAS1", "CAUSA_MULT"]

    return columns


def read_data(year, test=False):
    file_path = get_file_path(year)
    columns = get_columns(year)

    if test:
        print_log_message("THIS IS A TEST - READING ONLY 1000 ROWS FOR {}".format(year))
        if year < 2018:
            reader = pd.read_stata(
                file_path, iterator=True, convert_categoricals=False, columns=columns
            )
            df = reader.get_chunk(1000)
        else:
            df = pd.read_csv(file_path, nrow=1000, usecols=columns)
    else:
        print_log_message("Reading data for {}".format(year))
        if year < 2018:
            df = pd.read_stata(file_path, convert_categoricals=False, columns=columns)
        else:
            df = pd.read_csv(file_path, encoding="latin-1", usecols=columns)

    df.columns = df.columns.str.lower()

    return df


def get_age_map(year):

    if year in range(1998, 2008):
        age_map = {
            "01": "Early Neonatal",
            "02": "Early Neonatal",
            "03": "Late Neonatal",
            "04": "1-5 months",
            "05": "1-5 months",
            "06": "6-11 months",
            "07": "12 to 23 months",
            "08": "2 to 4",
            "09": "5 to 9",
            "10": "10 to 14",
            "11": "15 to 19",
            "12": "20 to 24",
            "13": "25 to 29",
            "14": "30 to 34",
            "15": "35 to 39",
            "16": "40 to 44",
            "17": "45 to 49",
            "18": "50 to 54",
            "19": "55 to 59",
            "20": "60 to 64",
            "21": "65 to 69",
            "22": "70 to 74",
            "23": "75 to 79",
            "24": "80 to 84",
            "25": "85 plus",
            "26": "Unknown",
        }
    elif year in range(2008, 2018):
        age_map = {
            "00": "Early Neonatal",
            "01": "Early Neonatal",
            "02": "Early Neonatal",
            "03": "Late Neonatal",
            "04": "1-5 months",
            "05": "1-5 months",
            "06": "6-11 months",
            "07": "12 to 23 months",
            "08": "2 to 4",
            "09": "5 to 9",
            "10": "10 to 14",
            "11": "15 to 19",
            "12": "20 to 24",
            "13": "25 to 29",
            "14": "30 to 34",
            "15": "35 to 39",
            "16": "40 to 44",
            "17": "45 to 49",
            "18": "50 to 54",
            "19": "55 to 59",
            "20": "60 to 64",
            "21": "65 to 69",
            "22": "70 to 74",
            "23": "75 to 79",
            "24": "80 to 84",
            "25": "85 to 89",
            "26": "90 to 94",
            "27": "95 plus",
            "28": "95 plus",
            "29": "Unknown",
        }
    elif year >= 2018:
        age_map = {
            0: "Early Neonatal",
            1: "Early Neonatal",
            2: "Early Neonatal",
            3: "Late Neonatal",
            4: "1-5 months",
            5: "1-5 months",
            6: "6-11 months",
            7: "12 to 23 months",
            8: "2 to 4",
            9: "5 to 9",
            10: "10 to 14",
            11: "15 to 19",
            12: "20 to 24",
            13: "25 to 29",
            14: "30 to 34",
            15: "35 to 39",
            16: "40 to 44",
            17: "45 to 49",
            18: "50 to 54",
            19: "55 to 59",
            20: "60 to 64",
            21: "65 to 69",
            22: "70 to 74",
            23: "75 to 79",
            24: "80 to 84",
            25: "85 to 89",
            26: "90 to 94",
            27: "95 plus",
            28: "95 plus",
            29: "Unknown",
        }

    return age_map


def create_age_group_id(df, year):

    age_map = get_age_map(year)

    gbd_ages = get_ages()
    gbd_ages = gbd_ages[["age_group_name", "age_group_id"]]
    assert ~gbd_ages.duplicated().any()

    df["age_group_name"] = df.gru_ed1.map(age_map)
    assert df.age_group_name.notnull().all()

    df = df.merge(gbd_ages, how="left", on="age_group_name")
    report_if_merge_fail(df, "age_group_id", "age_group_name")
    return df


def create_sex_id(df):

    sex_map = {1: 1, 2: 2, 3: 9}

    df["sex_id"] = df.sexo.astype(int).map(sex_map)
    assert df.sex_id.notnull().all()
    return df


def create_nid(df, year):

    year_to_nid = {
        1998: "NUMBER",
        1999: "NUMBER",
        2000: "NUMBER",
        2001: "NUMBER",
        2002: "NUMBER",
        2003: "NUMBER",
        2004: "NUMBER",
        2005: "NUMBER",
        2006: "NUMBER",
        2007: "NUMBER",
        2008: "NUMBER",
        2009: "NUMBER",
        2010: "NUMBER",
        2011: "NUMBER",
        2012: "NUMBER",
        2013: "NUMBER",
        2014: "NUMBER",
        2015: "NUMBER",
        2016: "NUMBER",
        2017: "NUMBER",
        2018: "NUMBER",
        2019: "NUMBER",
        2020: "NUMBER",
        2021: "NUMBER",
        2022: "NUMBER",
    }

    df["nid"] = year_to_nid[year]

    return df


def mark_hospital_deaths(df):
    df.sit_defun = df.sit_defun.astype(int)
    df["hospital_death"] = 0
    df.loc[df.sit_defun == 1, "hospital_death"] = 1
    return df


def clean_cause_column(df, column):
    not_like_ICD = ~df[column].str.match("(^[A-Z][0-9]{2,4}$)|(^0000$)")
    df.loc[not_like_ICD, column] = df[column].apply(lambda x: x[:-1])

    not_like_ICD = ~df[column].str.match("(^[A-Z][0-9]{2,4}$)|(^0000$)")
    not_like_ICD_list = df.loc[not_like_ICD, column].tolist()
    indices_changed = df.loc[not_like_ICD].index.tolist()
    df.loc[not_like_ICD, column] = "0000"

    return df, not_like_ICD_list, indices_changed


def create_underlying_cause(df, year, verbose=False):
    if year < 2018:
        df["cause"] = df.c_bas1.str.rstrip("X").str.upper()

        (df, ws, ic) = clean_cause_column(df, "cause")

        if verbose:
            print_log_message(
                "Found the following unexpected causes in cause column: {}".format(set(ws))
            )
            print_log_message("Changed {} deaths".format(len(set(ic))))
    else:
        df["cause"] = df["c_bas1"]

    validate_cause_column(df, "cause")

    missing_ucod = df["cause"] == "0000"
    assert (
        len(df[missing_ucod]) < 5
    ), "There are {} rows missing an underlying cause of death!".format(len(df[missing_ucod]))
    df = df[~missing_ucod]

    return df


def create_multiple_causes(df, year, verbose=False):
    if year < 2019:
        if year in range(1998, 2008):
            cols = ["c_dir1", "c_ant1", "c_ant2", "c_ant3", "c_pat1"]
        elif year in range(2008, 2019):
            cols = [
                "c_dir1",
                "c_dir12",
                "c_ant1",
                "c_ant12",
                "c_ant2",
                "c_ant22",
                "c_ant3",
                "c_ant32",
                "c_pat1",
                "c_pat2",
            ]

        unexpected = []
        indices_changed = []
        for col in cols:
            suffix = col.split("_")[1]
            mc_col = "multiple_cause_" + suffix
            df[mc_col] = df[col].str.rstrip("X").str.upper()
            df.loc[df[mc_col] == "", mc_col] = "0000"

            (df, ws, ic) = clean_cause_column(df, mc_col)
            unexpected += ws
            indices_changed += ic
            validate_cause_column(df, mc_col)

        if verbose:
            print_log_message(
                "Found the following unexpected causes in the multiple_cause columns: {}".format(
                    set(unexpected)
                )
            )
            print_log_message("Changed {} deaths".format(len(set(indices_changed))))
    else:
        df["causa_mult"] = df["causa_mult"].str.replace(" ", "/")
        df["causa_mult"] = df["causa_mult"].str.replace("\\", "/")
        df["causa_mult"] = df["causa_mult"].str.replace("*", "/")

        df["counts"] = df["causa_mult"].str.count("/")
        num_split = df["counts"].max()

        mc_cols = ["multiple_cause_" + str(i) for i in range(1, num_split + 2)]
        df[mc_cols] = df["causa_mult"].str.split("/", n=num_split, expand=True)

        for col in mc_cols:
            assert df[col].str.count("/").any() == 0

        for col in mc_cols:
            df.loc[df[col].isna(), col] = "0000"
            df.loc[df[col] == "None", col] = "0000"
            df[col] = df[col].str.upper()
            df.loc[~df[col].str.len().isin(range(3, 7)), col] = df[col].str[0:4]
            df.loc[~df[col].str.match("(^[A-Z][0-9]{2,4}$)|(^0000$)"), col] = df[col].str[:-1]
            bad = df.loc[~df[col].str.match("(^[A-Z][0-9]{2,4}$)|(^0000$)"), col].unique().tolist()
            if verbose:
                print(bad)
            assert len(bad) < 10
            df.loc[df[col].str.match("^([A-Z][0-9]{2,}|0000)$") == False, col] = "0000"
            validate_cause_column(df, col)

    return df


def handle_part_II(df, drop_p2=None):
    assert drop_p2 is not None, "Need to know whether to keep part II"

    p2_cols = [col for col in df.columns if col in ["multiple_cause_pat1", "multiple_cause_pat2"]]

    if drop_p2:
        df.drop(p2_cols, axis="columns", inplace=True)
    else:
        mc_cols = [col for col in df.columns if "multiple_cause" in col]
        for mc_col in mc_cols:
            suffix = mc_col.split("_")[-1]
            indicator_col = "pII_" + suffix
            if mc_col in p2_cols:
                df[indicator_col] = 1
            else:
                df[indicator_col] = 0
    return df


def clean_col_dane(year, drop_p2=None, test=False):
    assert year in YEARS, "No data available for {}".format(year)

    df = read_data(year, test=test)

    print_log_message("Cleaning data for {}".format(year))

    df = create_age_group_id(df, year)
    df = create_sex_id(df)
    df["year_id"] = year
    df["location_id"] = 125
    df["code_system_id"] = 1
    df = create_nid(df, year)
    df = mark_hospital_deaths(df)

    df = create_underlying_cause(df, year)
    df = create_multiple_causes(df, year)

    df = handle_part_II(df, drop_p2=drop_p2)

    df["deaths"] = 1

    mcodformatter = MultipleCauseFormatter(df, source=SOURCE, data_type_id=9, drop_p2=drop_p2)
    df = mcodformatter.get_formatted_dataframe()

    return df


if __name__ == "__main__":
    test = False
    write_metadata = False
    drop_p2 = True
    add_to_project = False
    project_ids = [1, 2, 4]

    for year in YEARS:
        df = clean_col_dane(year, drop_p2=drop_p2, test=test)
        if write_metadata:
            mcodformatter = MultipleCauseFormatter(
                df, source=SOURCE, data_type_id=9, drop_p2=drop_p2
            )
            mcodformatter.write_metadata(df)

        if add_to_project:
            mcodformatter = MultipleCauseFormatter(
                df, source=SOURCE, data_type_id=9, value_cols=VALUE_COLS, drop_p2=True
            )
            assert type(project_ids) == list, "Error: project_ids need to be a list"
            mcodformatter.add_to_project(df, project_ids)
