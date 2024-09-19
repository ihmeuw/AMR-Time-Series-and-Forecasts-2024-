import re
from builtins import zip

import numpy as np
import pandas as pd

from cod_prep.claude.configurator import Configurator
from cod_prep.downloaders import (
    add_age_metadata,
    add_cause_metadata,
    add_code_metadata,
    get_cause_map,
    get_cause_package_map,
    get_current_cause_hierarchy,
)
from cod_prep.utils import (
    clean_icd_codes,
    print_log_message,
    report_duplicates,
    report_if_merge_fail,
)
from mcod_prep.utils.causes import (
    get_all_related_syndromes,
    get_child_to_available_parent_syndrome,
    get_infsyn_hierarchy,
)
from mcod_prep.utils.misc import clean_col_strings


class MCoDMapper:

    cache_options = {"force_rerun": True, "block_rerun": False}
    conf = Configurator()
    inj_causes = ["x59", "y34"]
    int_cause_name_dict = {
        "pulmonary_embolism": ["pulmonary embolism"],
        "right_hf": ["right heart failure"],
        "left_hf": ["left heart failure"],
        "unsp_hf": ["unspecified heart failure"],
        "arterial_embolism": ["arterial embolism"],
        "aki": ["acute kidney failure"],
        "x59": ["unspecified external factor x59"],
        "y34": ["external causes udi,type unspecified-y34"],
        "drug_overdose": ["drug_overdose"],
        "hepatic_failure": ["hepatic failure"],
        "empyema": ["pleurisy,pyothorax"],
        "arf": ["acute respiratory failure"],
        "pneumonitis": "",
        "unsp_cns": ["unspecified cns signs and symptom"],
        "gi_bleeding": ["gastrointestinal bleeding"],
        "arrhythmia": ["cardiac rhythm disorders"],
        "amyloidosis": "",
        "peritonitis": ["peritonitis&acute abdomen"],
        "atherosclerosis": "",
        "fe_acid_base": ["fluid, electrolit,acid base disorders"],
        "hypertension": "",
        "crf": ["chronic respiratory failure"],
        "pneumothorax": "",
        "cerebral_palsy": ["cerebral palsy"],
        "plegia": "",
        "cachexia": "",
        "cardiac_arrest": ["shock,cardiac arrest,coma"],
        "osteomyelitis": "",
        "alc_hepatic_failure": ["alcoholic hepatic failure"],
        "sepsis": "",
        "renal_failure": ["renal failure"],
        "pulmonary_embolism_non_septic": "non septic pulmonary embolism",
        "hf": ["right heart failure", "left heart failure", "unspecified heart failure"],
        "embolism": ["arterial embolism", "pulmonary embolism"],
    }
    infectious_syndromes = [
        "infectious_syndrome",
        "pathogen",
        "L2_meningitis",
        "L2_encephalitis",
        "L2_myelitis_meningoencephalitis_other",
        "L2_skin_infection",
        "L2_oral_infection",
        "L2_eye_infection",
        "L2_diarrhea",
        "L2_typhoid_paratyphoid_ints",
        "L2_hepatitis",
        "L2_gastrointestinal_infection_other",
        "L2_genital_infection",
        "L2_sexually_transmitted_infection",
        "L2_blood_stream_infection",
        "L2_endocarditis",
        "L2_tb",
        "L2_leprosy",
        "L2_other_parasitic_infection",
        "L2_unspecified_site_infection",
        "L1_bone_joint_infection",
        "L1.5_myocarditis_pericarditis_carditis",
        "L1_peritoneal_and_intra_abdomen_infection",
        "L2_lower_respiratory_infection_hosp",
        "L2_lower_respiratory_infection_comm",
        "L2_upper_respiratory_infection",
        "L2_mix_respiratory_infection",
        "L2_urinary_tract_infection_comm",
        "L2_urinary_tract_infection_hosp",
        "L2_OAMOIS",
    ]

    possible_int_causes = (
        list(int_cause_name_dict.keys())
        + infectious_syndromes
        + gbd_acauses
        + ["custom", "cross_tabulate", "standard_map", "flu_type"]
    )
    combination_syndrome_dict = {
        "L2_myocarditis": "L1.5_myocarditis_pericarditis_carditis",
        "L2_pericarditis": "L1.5_myocarditis_pericarditis_carditis",
        "L2_carditis": "L1.5_myocarditis_pericarditis_carditis",
        "L2_bone_joint_bacterial_inf": "L1_bone_joint_infection",
        "L2_bone_joint_non_bacterial_inf": "L1_bone_joint_infection",
        "L2_peritonitis": "L1_peritoneal_and_intra_abdomen_infection",
        "L2_intra_abdomen_abscess": "L1_peritoneal_and_intra_abdomen_infection",
    }
    remainder_syndromes = [
        "L1_non_infectious_cause",
        "L2_inj_medical",
        "L2_other_injuries",
        "L2_non_infectious_non_inj",
        "L2_non_specific_icd_code",
    ]
    null_code = "0000"

    def __init__(
        self,
        int_cause,
        drop_p2,
        code_system_id=None,
        code_map_version_id=None,
        path_to_ucause_map=None,
        path_to_int_cause_map=None,
        infsyn_set_version="gram2",
    ):
        self.int_cause = int_cause
        self.drop_p2 = drop_p2

        self.code_system_id = code_system_id
        self.code_map_version_id = code_map_version_id
        self.path_to_ucause_map = path_to_ucause_map
        self.path_to_int_cause_map = path_to_int_cause_map

        self.infsyn_set_version = infsyn_set_version

        if self.code_system_id and self.code_map_version_id:
            self.is_real_code_system = True
        elif self.path_to_ucause_map or self.path_to_int_cause_map:
            assert self.int_cause == "infectious_syndrome", (
                "You specified paths to custom maps but you are not trying"
                "to capture infectious syndromes - this option is not implemented"
            )
            self.is_real_code_system = False
        else:
            raise AssertionError(
                "You must pass either code system and code map version "
                "or paths to ucause and int_cause maps"
            )

        if self.infsyn_set_version not in ["v1", "gram2"]:
            assert self.int_cause == "infectious_syndrome", (
                "You specified a non-default infectious syndrome hierarchy "
                "but you are not trying to capture infectious syndromes - "
                "this option is not implemented"
            )
            assert not self.is_real_code_system, (
                "You specified a non-default infectious syndrome hiearchy "
                "but you are attempting to use code-system based maps - "
                "ICD codes have not been mapped to non-standard hierarchies"
            )

        assert (
            self.int_cause in self.possible_int_causes
        ), f"{self.int_cause} is not a valid intermediate cause"
        try:
            self.full_cause_name = self.int_cause_name_dict[self.int_cause]
        except KeyError:
            self.full_cause_name = self.int_cause
        if self.full_cause_name == "":
            self.full_cause_name = self.int_cause
        if type(self.full_cause_name) != list:
            self.full_cause_name = [self.full_cause_name]

    @staticmethod
    def get_code_columns(df):
        col_names = list(df.columns)
        code_cols = [x for x in col_names if "multiple_cause" in x and "pII" not in x]
        if "cause" in df:
            code_cols += ["cause"]
        return code_cols

    @staticmethod
    def _get_cause_num(mcod_col):
        if mcod_col.startswith("cause"):
            return "0"
        else:
            assert re.match(
                r"^multiple_cause_[a-z]*[0-9]*", mcod_col
            ), f"column {mcod_col} does not match expected format: multiple_cause_x"
            return mcod_col.split("_")[2]

    @staticmethod
    def prep_raw_mapped_cause_dictionary(raw_cols, mapped_cols):
        raw_cols = sorted(raw_cols, key=MCoDMapper._get_cause_num)
        mapped_cols = sorted(mapped_cols, key=MCoDMapper._get_cause_num)
        return dict(list(zip(raw_cols, mapped_cols)))

    @staticmethod
    def fix_icd_codes(df, codes, code_system_id):
        df[codes] = df[codes].astype(str)
        if code_system_id == 6:
            df.loc[df["cause"].str.contains("^[89]"), "cause"] = "E" + df["cause"]
        return df

    @staticmethod
    def prep_cause_package_map(cause_package_map):
        check_map = cause_package_map[["map_id", "map_type"]].drop_duplicates()
        report_duplicates(check_map, "map_id")
        return cause_package_map.set_index("value")["map_id"].to_dict()

    @staticmethod
    def prep_crosstab_cause_package_map(cause_package_map):
        cpm = cause_package_map[["value", "map_id", "map_type"]].drop_duplicates()
        report_duplicates(cpm, "value")
        assert cpm.map_type.isin(["cause_id", "package_id"]).all()
        map_type_to_prefix = {"cause_id": "c_", "package_id": "p_"}
        cpm["map_id"] = cpm.apply(
            lambda x: map_type_to_prefix[x["map_type"]] + str(int(x["map_id"])), axis=1
        )
        return cpm.set_index("value")["map_id"].to_dict()

    @staticmethod
    def prep_cause_map(cause_map):
        cause_map["value"] = clean_icd_codes(cause_map["value"], remove_decimal=True)
        cause_map = cause_map.drop_duplicates(["code_system_id", "value"])
        cause_map["code_id"] = cause_map["code_id"].astype(int)
        return cause_map.set_index("value")["code_id"].to_dict()

    def map_cause_codes(self, df, coi_map, coi, cols_to_map=None):
        if not cols_to_map:
            cols_to_map = MCoDMapper.get_code_columns(df)

        for col in cols_to_map:
            df[col] = df[col].fillna(MCoDMapper.null_code)
            df[col] = df[col].astype(object)
            if self.code_system_id not in [
                1,
                6,
            ]:
                df = clean_col_strings(df, col, deep_clean=False)
            df[col + "_" + coi] = df[col].map(coi_map)
            df.loc[df[col] == MCoDMapper.null_code, col + "_" + coi] = MCoDMapper.null_code

        return df

    @staticmethod
    def trim_and_remap(self, df, code_dict, cause_map):
        df = df.copy()

        for n in reversed(range(3, 6)):
            for code, mapped_code in list(code_dict.items()):
                temp_code = "temp_" + code
                df[temp_code] = df[code].copy()
                try:
                    df.loc[df[mapped_code].isnull(), temp_code] = df[temp_code].apply(
                        lambda x: x[0:n]
                    )
                except TypeError:
                    if mapped_code != "cause_mapped":
                        df[mapped_code] = MCoDMapper.null_code
                    else:
                        print("problem code here..." + df[code])
                df.loc[df[mapped_code].isnull(), mapped_code] = df[temp_code].map(cause_map)
                df = df.drop(temp_code, axis=1)
        return df

    def fill_missing_ICDs(self, coi_map, map_type):
        if self.code_system_id == 1:
            if map_type == "b_cause_id":
                coi_map.update(
                    {
                        "acause_inj_trans_road_4wheel": 693,
                        "H63": 2136,
                        "T896": 2136,
                        "T771": 2136,
                        "Y67": 919,
                        "C89": 919,
                        "0000": 919,
                        "Y6899": 919,
                        "Y39199": 919,
                        "none": 919,
                        "K15": 2136,
                        "K150": 2136,
                        "K151": 2136,
                        "K153": 2136,
                        "T8901": 2136,
                        "T8903": 2136,
                        "T8900": 2136,
                        "T8902": 2136,
                        "T775": 2136,
                        "T725": 2136,
                        "T72": 2136,
                        "T77": 2136,
                        "H64": 1061,
                        "I91": 507,
                        "I911": 507,
                        "T720": 2136,
                        "T721": 2136,
                        "T722": 2136,
                        "T728": 2136,
                        "T729": 2136,
                        "T779": 2136,
                        "T89": 2136,
                        "T892": 2136,
                        "T897": 2136,
                        "T894": 2136,
                        "T899": 2136,
                        "T99": 2136,
                        "T994": 2136,
                        "T999": 2136,
                        "T896": 2136,
                        "Q570": 919,
                        "Q579": 919,
                        "Q589": 2136,
                        "Q585": 2136,
                        "Q588": 2136,
                        "Q592": 2136,
                    }
                )
            elif map_type == "code_id":
                coi_map.update(
                    {
                        "acause_inj_trans_road_4wheel": 35176,
                        "acause_inj_fires": 103630,
                        "H63": 38783,
                        "T896": 38783,
                        "T771": 38783,
                        "Y67": 95265,
                        "C89": 95265,
                        "0000": 95265,
                        "Y6899": 95265,
                        "Y39199": 95265,
                        "none": 95265,
                        "K15": 38783,
                        "K150": 38783,
                        "K151": 38783,
                        "K153": 38783,
                        "T8901": 38783,
                        "T8903": 38783,
                        "T8900": 38783,
                        "T8902": 38783,
                        "T775": 38783,
                        "T725": 38783,
                        "T72": 38783,
                        "T77": 38783,
                        "H64": 26902,
                        "I91": 11123,
                        "I911": 11123,
                        "T720": 38783,
                        "T721": 38783,
                        "T722": 38783,
                        "T728": 38783,
                        "T729": 38783,
                        "T779": 38783,
                        "T89": 38783,
                        "T892": 38783,
                        "T897": 38783,
                        "T894": 38783,
                        "T899": 38783,
                        "T99": 38783,
                        "T994": 38783,
                        "T999": 38783,
                        "T896": 38783,
                        "Q570": 95265,
                        "Q579": 95265,
                        "Q589": 38783,
                        "Q585": 38783,
                        "Q588": 38783,
                        "Q592": 38783,
                    }
                )
            elif map_type == "infectious_syndrome":
                coi_map.update(
                    {
                        "acause_inj_trans_road_4wheel": "L2_other_injuries",
                        "H63": "L2_non_specific_icd_code",
                        "T896": "L2_non_specific_icd_code",
                        "T771": "L2_non_specific_icd_code",
                        "Y6899": "0000",
                        "Y39199": "0000",
                        "Y67": "0000",
                        "C89": "0000",
                        "K15": "L2_non_infectious_non_inj",
                        "K150": "L2_non_infectious_non_inj",
                        "K151": "L2_non_infectious_non_inj",
                        "K15": "L2_non_infectious_non_inj",
                        "K153": "L2_non_infectious_non_inj",
                        "T775": "L2_non_specific_icd_code",
                        "T725": "L2_non_specific_icd_code",
                        "T72": "L2_non_specific_icd_code",
                        "T77": "L2_non_specific_icd_code",
                        "T779": "L2_non_specific_icd_code",
                        "T771": "L2_non_specific_icd_code",
                        "acause_inj_fires": "L2_other_injuries",
                        "T89": "L2_non_specific_icd_code",
                        "T8901": "L2_non_specific_icd_code",
                        "T8903": "L2_non_specific_icd_code",
                        "T899": "L2_non_specific_icd_code",
                        "T8903": "L2_non_specific_icd_code",
                        "T894": "L2_non_specific_icd_code",
                        "T8902": "L2_non_specific_icd_code",
                        "Q589": "L2_non_infectious_non_inj",
                        "Q588": "L2_non_infectious_non_inj",
                        "Q592": "L2_non_infectious_non_inj",
                        "Q585": "L2_non_infectious_non_inj",
                        "none": "0000",
                    }
                )
            else:
                print("Not a valid map_type. Please use 'b_cause' or 'code_id'.")
        elif self.code_system_id == 6:
            if map_type == "b_cause_id":
                coi_map.update(
                    {
                        "V3812": 693,
                        "V942": 695,
                        "0000": 919,
                        "none": 919,
                    }
                )
            elif map_type == "code_id":
                coi_map.update(
                    {
                        "V3812": 64155,
                        "V942": 63998,
                        "0000": 95273,
                        "none": 95273,
                    }
                )
            elif map_type == "infectious_syndrome":
                coi_map.update(
                    {
                        "none": "0000",
                        "V3812": "L2_non_specific_icd_code",
                        "V942": "L2_non_specific_icd_code",
                    }
                )
            else:
                print("Not a valid map_type. Please use 'b_cause_id' or 'code_id'.")
        else:
            print("Not a valid code system")
        return coi_map

    @staticmethod
    def assert_all_mapped(df, code_dict):
        unmapped_dfs = []
        for code, mapped_code in code_dict.items():
            unmapped = df.loc[df[mapped_code].isnull(), [code]].drop_duplicates()
            unmapped = unmapped.rename(columns={code: "code_value"})
            unmapped["column"] = code
            unmapped_dfs.append(unmapped)
        unmapped = pd.concat(unmapped_dfs)
        if len(unmapped) > 0:
            raise AssertionError(f"The following codes failed to be mapped: \n{unmapped}")

    def drop_duplicated_chain(df, chain_cols, suffix=None):
        if suffix != None:
            ucod_col = f"cause_{suffix}"
        else:
            ucod_col = "cause"

        if ucod_col in chain_cols:
            chain_cols = [x for x in chain_cols if x != ucod_col]

        for col in chain_cols:
            df.loc[df[col] == df[ucod_col], col] = MCoDMapper.null_code

            check_chain = [x for x in chain_cols if x != col]
            for check_col in check_chain:
                df.loc[
                    (df[check_col] != MCoDMapper.null_code) & (df[check_col] == df[col]), check_col
                ] = MCoDMapper.null_code

        return df

    def prep_int_cause_map(self):
        code_system_name = {1: "icd10", 6: "icd9"}.get(self.code_system_id, None)
        if self.int_cause == "sepsis":
            df = pd.read_csv(
                "FILEPATH",
                dtype={"icd_code": object},
            )
            df["icd_code"] = clean_icd_codes(df["icd_code"], remove_decimal=True)
            df = df[["icd_code", "sepsis"]].drop_duplicates()
            report_duplicates(df, "icd_code")
            mcod_map = dict(list(zip(df["icd_code"], df["sepsis"])))
            assert set(mcod_map.values()).issubset(
                {
                    "ucod_no need of-implicit",
                    "explicit",
                    "ucod_need of-implicit",
                    "no_sepsis",
                    "implicit-of-code",
                }
            )

        elif self.int_cause in self.infectious_syndromes:
            if self.is_real_code_system:
                df = pd.read_csv(
                    "FILEPATH",
                    dtype={"icd_code": object},
                )
                if self.code_system_id in [1, 6]:
                    df["icd_code"] = clean_icd_codes(df["icd_code"], remove_decimal=True)
                df = df[["icd_code", "infectious_syndrome"]].drop_duplicates()
                mcod_map = dict(list(zip(df["icd_code"], df["infectious_syndrome"])))
            else:
                df = pd.read_csv("FILEPATH", encoding="latin1")[["cause", "infectious_syndrome"]]
                df = clean_col_strings(df, "cause", deep_clean=False)
                df = df.drop_duplicates()
                report_duplicates(df, "cause")
                mcod_map = dict(list(zip(df["cause"], df["infectious_syndrome"])))

        elif self.int_cause == "custom":
            df = pd.read_excel(
                "FILEPATH",
                dtype={"icd_code": str},
            )
            df["icd_code"] = clean_icd_codes(df["icd_code"], remove_decimal=True)
            report_duplicates(df, ["icd_code", "code_system"])
            df = df.query(f"code_system == '{code_system_name}'")
            assert len(df) > 0, "Your custom map has no valid mappings"
            df["map_col"] = "custom"
            mcod_map = dict(list(zip(df["icd_code"], df["map_col"])))

        elif self.int_cause in self.gbd_acauses:
            df = get_cause_map(code_map_version_id=self.code_map_version_id, **self.cache_options)
            df = df.loc[df.value != "E90.0"]
            df = df.loc[df.value.str[-1] != "."]
            df["icd_code"] = clean_icd_codes(df["value"], remove_decimal=True)
            df = add_cause_metadata(
                df,
                "acause",
                cause_set_version_id=self.conf.get_id("cause_set_version"),
                **self.cache_options,
            )
            report_if_merge_fail(df, "acause", "cause_id")
            report_duplicates(df, ["icd_code"])
            assert len(df) > 0, "Your specified GBD acause has no valid mappings"
            mcod_map = df.set_index("icd_code")["acause"].to_dict()

        elif self.int_cause == "cross_tabulate" or self.int_cause == "standard_map":

            hospital_claims_map = pd.read_csv(
                "FILEPATH",
                encoding="utf-8-sig",
            )
            bcause_map = hospital_claims_map.loc[
                hospital_claims_map["code_system"] == self.code_system_id
            ]

            bcause_map = bcause_map[["icd_code", "b_cause_id"]]
            bcause_map["icd_code"] = clean_icd_codes(bcause_map["icd_code"], remove_decimal=True)

            bcause_map["b_cause_id"] = bcause_map["b_cause_id"].astype(int)

            mcod_map = bcause_map.set_index("icd_code")["b_cause_id"].to_dict()
            return mcod_map

        elif self.int_cause == "flu_type":
            df = pd.read_csv("FILEPATH", dtype={"icd_code": str})
            df = df.query(f"code_system_id == {self.code_system_id}")
            assert not df.empty
            df["icd_code"] = clean_icd_codes(df["icd_code"], remove_decimal=True)
            mcod_map = df.set_index("icd_code", verify_integrity=True)["flu_type"].to_dict()

        else:
            df = pd.read_csv("FILEPATH", encoding="latin-1")
            df = df[["icd_code", "package_description", "code_system"]].drop_duplicates()
            df["icd_code"] = clean_icd_codes(df["icd_code"], remove_decimal=True)
            df["package_description"] = df["package_description"].astype(str).str.lower()
            df["code_system"] = df["code_system"].str.lower()
            df = df.loc[df["package_description"].isin(self.full_cause_name)]
            report_duplicates(df, ["icd_code", "code_system"])
            df = df.query(f'code_system == "{code_system_name}"')
            assert (
                len(df) > 0
            ), f"There are no mappings for {code_system_name}, {self.full_cause_name}"
            mcod_map = dict(list(zip(df["icd_code"], df["package_description"])))

        return mcod_map

    def capture_int_cause(self, df, int_cause_cols):
        if self.int_cause == "sepsis":
            df[int_cause_cols] = df[int_cause_cols].fillna("no_sepsis")
            df.loc[(df[int_cause_cols] == "explicit").any(axis=1), "sepsis"] = "explicit"
            not_explicit = df.sepsis != "explicit"
            ucod_no_need_of = df.cause_sepsis == "ucod_no need of-implicit"
            ucod_need_of = df.cause_sepsis == "ucod_need of-implicit"
            has_of = (df[int_cause_cols] == "implicit-of-code").any(axis=1)
            fatal = df.deaths > 0
            df.loc[
                not_explicit & (ucod_no_need_of | (ucod_need_of & has_of) | (ucod_need_of & fatal)),
                "sepsis",
            ] = "implicit"

            df["sepsis"] = df["sepsis"].fillna("no_sepsis")
        elif self.int_cause in self.inj_causes:
            df = self.capture_injuries_pattern(df, int_cause_cols)
        elif self.int_cause in self.infectious_syndromes:
            df = self.apply_syndrome_severity(df, int_cause_cols, self.int_cause)
            if self.int_cause == "pathogen":
                df = self.extract_pathogens(df, int_cause_cols, self.int_cause)
        elif self.int_cause == "cross_tabulate":
            df = self.cross_tabulate_shaping(df, int_cause_cols)
        elif self.int_cause == "standard_map":
            df = self.standard_map_prep(df, int_cause_cols)
        elif self.int_cause == "flu_type":
            df[int_cause_cols] = df[int_cause_cols].fillna("non_flu")
            df[self.int_cause] = df[int_cause_cols].apply(
                lambda x: ",".join(sorted(set(x) - {"non_flu", "0000"})), axis=1
            )
            df[self.int_cause] = df[self.int_cause].replace("", "non_flu")
        else:
            for col in int_cause_cols:
                df[col] = df[col].fillna("other")
                df.loc[df[col].isin(self.full_cause_name), self.int_cause] = 1
            df[self.int_cause] = df[self.int_cause].fillna(0)

        print_log_message("Finished capturing int cause")
        assert df[self.int_cause].notnull().values.all()
        return df

    def extract_pathogens(self, df, int_cause_cols: list, int_cause: str):
        raw_cause_cols = [
            x for x in df.columns if "cause" in x and "pathogen" not in x and x != "cause_id"
        ]

        infsyn = get_infsyn_hierarchy(infsyn_set_version=self.infsyn_set_version)

        infsyn = infsyn.drop_duplicates(subset=["infectious_syndrome"], keep="first")

        path_map = infsyn.loc[infsyn["level"] == 3, ["infectious_syndrome", "pathogen"]]
        path_map["pathogen"] = path_map["pathogen"].fillna(MCoDMapper.null_code)
        path_map = dict(list(zip(path_map["infectious_syndrome"], path_map["pathogen"])))

        df = MCoDMapper.map_cause_codes(self, df, path_map, "etiology", int_cause_cols)

        for col in [x for x in df.columns if "etiology" in x]:
            df[col] = df[col].fillna(MCoDMapper.null_code)
            df.rename(columns={col: col.replace("_pathogen", "")}, inplace=True)
        path_int_cause_cols = [x for x in df.columns if "etiology" in x]

        int_cause_col_dict = MCoDMapper.prep_raw_mapped_cause_dictionary(
            raw_cause_cols, path_int_cause_cols
        )

        print_log_message("Trimming ICD codes and remapping chain causes")
        df = MCoDMapper.trim_and_remap(df, int_cause_col_dict, path_map)

        print_log_message("Dropping duplicates after mapping for etiologies")
        eti_cols = [x for x in df.columns if "multiple_cause_" in x and "etiology" in x]
        df = MCoDMapper.drop_duplicated_chain(df, eti_cols, "etiology")

        cause_map = get_cause_map(
            code_map_version_id=self.code_map_version_id, **self.cache_options
        )
        cause_map["value"] = clean_icd_codes(cause_map["value"], remove_decimal=True)
        cause_map = cause_map.drop_duplicates(["code_system_id", "value"])
        cause_map = dict(list(zip(cause_map["value"], cause_map["cause_id"])))

        df = MCoDMapper.map_cause_codes(self, df, cause_map, "disease", raw_cause_cols)

        for col in [x for x in df.columns if "disease" in x]:
            df[col] = df[col].fillna(MCoDMapper.null_code)
            df.rename(columns={col: col.replace("_pathogen", "")}, inplace=True)
        disease_int_cause_cols = [x for x in df.columns if "disease" in x]

        int_cause_col_dict = MCoDMapper.prep_raw_mapped_cause_dictionary(
            raw_cause_cols, disease_int_cause_cols
        )

        print_log_message("Trimming ICD codes and remapping chain causes")
        df = MCoDMapper.trim_and_remap(df, int_cause_col_dict, cause_map)

        print_log_message("Dropping duplicates after mapping for diseases")
        dis_cols = [x for x in df.columns if "multiple_cause_" in x and "disease" in x]
        df = MCoDMapper.drop_duplicated_chain(df, dis_cols, "disease")

        return df

    def apply_syndrome_severity(self, df, int_cause_cols: list, syndrome: str):
        assert {
            "_row_id",
            "pathogen_from_cause",
            "infectious_syndrome",
            "l2_syndrome",
        }.isdisjoint(set(df))
        infsyn = get_infsyn_hierarchy(infsyn_set_version=self.infsyn_set_version)

        infsyn = infsyn.drop_duplicates(subset=["infectious_syndrome"], keep="first")
        severity_map = infsyn.set_index("infectious_syndrome")["severity"].to_dict()
        has_ucause = "cause" in df
        if has_ucause:
            cause_col_dict = {
                MCoDMapper._get_cause_num(cause_col): cause_col for cause_col in int_cause_cols
            }
            ucod_col = cause_col_dict["0"]
        true_syndromes = list(set(severity_map.keys()) - set(MCoDMapper.remainder_syndromes))

        for col in int_cause_cols:
            assert (
                df[col]
                .fillna("L1_non_infectious_cause")
                .isin(set(severity_map.keys()).union({MCoDMapper.null_code}))
                .all()
            )

        df["all_possible_syndromes"] = df[int_cause_cols].apply(
            lambda row: ",".join(
                set(
                    [
                        x
                        for x in sorted(row.values.astype(str))
                        if x
                        not in [
                            "0000",
                            "L1_non_infectious_cause",
                            "L2_non_specific_icd_code",
                            "L2_other_injuries",
                            "L2_non_infectious_non_inj",
                            "L2_inj_medical",
                            "none",
                        ]
                    ]
                )
            ),
            axis=1,
        )

        reverse_severity = {v: k for k, v in severity_map.items()}
        df = df.assign(
            infectious_syndrome=lambda x: x[int_cause_cols]
            .apply(lambda y: y.map(severity_map))
            .min(axis=1)
            .map(reverse_severity)
            .fillna("L1_non_infectious_cause")
        )

        assert df["infectious_syndrome"].notnull().all()

        bsi = get_all_related_syndromes("L2_blood_stream_infection", infsyn=infsyn)
        if has_ucause:
            df.loc[
                (df[ucod_col].isin(true_syndromes)) & (~df[ucod_col].isin(bsi)),
                "infectious_syndrome",
            ] = df[ucod_col]

        neonatal_bsi = get_all_related_syndromes("L3_bsi_neonatal_bacterial", infsyn=infsyn)
        df = add_age_metadata(df, ["age_group_days_end"], **self.cache_options)
        report_if_merge_fail(df, "age_group_days_end", "age_group_id")
        neonatal = df.age_group_days_end <= 27
        has_bsi = df[int_cause_cols].isin(bsi).any(axis=1)

        df.loc[(neonatal & has_bsi), "infectious_syndrome"] = (
            df[int_cause_cols]
            .apply(lambda x: list(set(x).intersection(bsi)), axis=1)
            .apply(
                lambda x: reverse_severity[min([severity_map[y] for y in x])] if len(x) > 0 else ""
            )
        )
        df = df.drop("age_group_days_end", axis="columns")

        syndrome_to_level_2_map = get_child_to_available_parent_syndrome(
            infsyn.query("level == 2").infectious_syndrome.unique(), infsyn=infsyn
        )
        level_2_to_children = {
            l2: get_all_related_syndromes(l2, infsyn=infsyn)
            for l2 in infsyn.query("level == 2").infectious_syndrome.unique()
        }
        df["l2_syndrome"] = df["infectious_syndrome"].map(syndrome_to_level_2_map)

        df["_row_id"] = list(range(0, df.shape[0]))
        syn_df = (
            df.loc[:, ["_row_id", "l2_syndrome"] + int_cause_cols]
            .melt(
                id_vars=["_row_id", "l2_syndrome"],
                value_vars=int_cause_cols,
                value_name="detailed_infsyn",
                var_name="chain_col",
            )
            .loc[:, ["_row_id", "l2_syndrome", "detailed_infsyn"]]
            .drop_duplicates()
        )
        syn_df["special_or_l1_syn"] = syn_df["l2_syndrome"].map(
            MCoDMapper.combination_syndrome_dict
        )
        syn_df.loc[syn_df["special_or_l1_syn"].notnull(), "l2_syndrome"] = syn_df[
            "special_or_l1_syn"
        ]
        syn_df.drop(columns=["special_or_l1_syn"], inplace=True)
        combination_syndrome_df = pd.DataFrame.from_dict(
            MCoDMapper.combination_syndrome_dict, orient="index"
        )
        combination_syndrome_df.reset_index(inplace=True)
        combination_syndrome_df.rename(columns={"index": "l2", 0: "new_combo"}, inplace=True)
        combination_syndrome_df["childrens_of_l2"] = combination_syndrome_df.apply(
            lambda row: get_all_related_syndromes(row["l2"], infsyn), axis=1
        )
        combination_syndrome_df = combination_syndrome_df.groupby("new_combo", as_index=False).agg(
            {"childrens_of_l2": "sum"}
        )
        cmb_syn_dict = dict(
            zip(combination_syndrome_df["new_combo"], combination_syndrome_df["childrens_of_l2"])
        )
        level_2_to_children.update(cmb_syn_dict)

        syn_df = pd.concat(
            [
                group.loc[group.detailed_infsyn.isin(level_2_to_children[key])]
                for key, group in syn_df.loc[
                    syn_df["l2_syndrome"].notnull(),
                ].groupby(["l2_syndrome"])
            ]
        )

        if self.int_cause == "infectious_syndrome":
            index_cols = list(set(df.columns).difference(int_cause_cols)) + ["cis"]
            df["cis"] = df["cause_infectious_syndrome"]
            value_cols = list(set(df.columns).intersection(["deaths", "admissions", "cases"]))
            orig_totals = df[value_cols].sum()
            df = pd.melt(
                df,
                id_vars=index_cols,
                value_vars=int_cause_cols,
                var_name="msyndrome_col",
                value_name="all_infectious_syndromes",
            )
            df = df.loc[
                df["all_infectious_syndromes"] != "0000",
            ]
            df["most_severe"] = 0
            df.loc[(df["all_infectious_syndromes"] == df["infectious_syndrome"]), "most_severe"] = 1
            df["infectious_syndrome"] = df["all_infectious_syndromes"]
            df["cause_infectious_syndrome"] = df["cis"]
            df.drop(columns=["all_infectious_syndromes", "msyndrome_col", "cis"], inplace=True)
            melt_orig_totals = df.loc[df["most_severe"] == 1, value_cols].sum()
            np.allclose(orig_totals, melt_orig_totals)

        if self.int_cause == "infectious_syndrome":
            syndrome_to_pathogen = infsyn.set_index("infectious_syndrome")["pathogen"].to_dict()
            syn_df["pathogen_from_cause"] = syn_df["detailed_infsyn"].map(syndrome_to_pathogen)

            sample_to_pathogens = (
                syn_df.loc[syn_df.pathogen_from_cause.notnull()]
                .groupby("_row_id")["pathogen_from_cause"]
                .aggregate(lambda x: ",".join(x.drop_duplicates().sort_values()))
                .to_dict()
            )

            dfsv = df.loc[
                df["most_severe"] == 1,
            ]
            dfsv["pathogen_from_cause"] = dfsv["_row_id"].map(sample_to_pathogens).fillna("none")
            dfnsv = df.loc[
                df["most_severe"] == 0,
            ]
            dfnsv = dfnsv.merge(
                infsyn[["infectious_syndrome", "pathogen"]], how="left", validate="many_to_one"
            )

            dfnsv.loc[dfnsv["pathogen"].isna(), "pathogen"] = "none"
            dfnsv.rename(columns={"pathogen": "pathogen_from_cause"}, inplace=True)
            dfnsv["l2_syndrome_dup"] = dfnsv["infectious_syndrome"].map(syndrome_to_level_2_map)
            dfnsv = dfnsv.loc[
                dfnsv["l2_syndrome"] != dfnsv["l2_syndrome_dup"],
            ]
            dfnsv.drop(columns=["l2_syndrome_dup"], inplace=True)
            df = pd.concat([dfsv, dfnsv])
            df = df.reset_index(drop=True)
            skin_infect_non_bacterial = (
                infsyn.loc[
                    infsyn["parent_infectious_syndrome"].isin(
                        ["L3_skin_viral", "L3_skin_fungal", "L3_skin_parasitic"]
                    ),
                    "infectious_syndrome",
                ]
                .unique()
                .tolist()
            )
            skin_infect_non_bacterial += ["L3_skin_viral", "L3_skin_fungal", "L3_skin_parasitic"]
            skin_pathogens_non_bacterial = (
                infsyn.loc[
                    infsyn["infectious_syndrome"].isin(skin_infect_non_bacterial), "pathogen"
                ]
                .unique()
                .tolist()
            )

            def subtract_lists_from_string(row, skin_pathogens_non_bacterial):
                row_list = row.split(",")

                reduced_list = [x for x in row_list if x not in skin_pathogens_non_bacterial]

                reduced_string = ",".join(reduced_list)
                return reduced_string

            if "deaths" in df.columns:
                filter_df = df.loc[
                    (df["deaths"] > 0) & (df["infectious_syndrome"].isin(skin_infect_non_bacterial))
                ]
                filter_df["pathogen_from_cause"] = filter_df["pathogen_from_cause"].apply(
                    lambda x: subtract_lists_from_string(x, skin_pathogens_non_bacterial)
                )
                df.update(filter_df[["pathogen_from_cause"]])
            df["pathogen_from_cause"] = df["pathogen_from_cause"].fillna("none")
            df.loc[df["pathogen_from_cause"] == "", "pathogen_from_cause"] = "none"

            df = df.loc[
                ~((df["most_severe"] == 0) & (df["pathogen_from_cause"] == "none")),
            ]

        else:
            if "_hosp" in self.int_cause:
                target_syndromes = get_all_related_syndromes(
                    self.int_cause.replace("_hosp", ""), infsyn
                )

            if "_comm" in self.int_cause:
                target_syndromes = get_all_related_syndromes(
                    self.int_cause.replace("_comm", ""), infsyn
                )

            if "_comm" not in self.int_cause and "_hosp" not in self.int_cause:
                if self.int_cause not in MCoDMapper.combination_syndrome_dict.values():
                    target_syndromes = get_all_related_syndromes(self.int_cause, infsyn)
                else:
                    target_syndromes = level_2_to_children[self.int_cause]

            samples_to_flag = syn_df.loc[
                syn_df.detailed_infsyn.isin(target_syndromes), "_row_id"
            ].unique()
            df[self.int_cause] = df._row_id.isin(samples_to_flag).astype(int)

            if "_comm" in self.int_cause:
                df.loc[
                    (df["infectious_syndrome"] != df[ucod_col]) & (df[self.int_cause] == 1),
                    self.int_cause,
                ] = 0

            if "_hosp" in self.int_cause:
                df.loc[
                    (df["infectious_syndrome"] == df[ucod_col]) & (df[self.int_cause] == 1),
                    self.int_cause,
                ] = 0

            if self.int_cause == "L2_skin_infection":
                skin_infect_non_bacterial = (
                    infsyn.loc[
                        infsyn["parent_infectious_syndrome"].isin(
                            ["L3_skin_viral", "L3_skin_fungal", "L3_skin_parasitic"]
                        ),
                        "infectious_syndrome",
                    ]
                    .unique()
                    .tolist()
                )
                skin_infect_non_bacterial += [
                    "L3_skin_viral",
                    "L3_skin_fungal",
                    "L3_skin_parasitic",
                ]
                skin_infect_bacterial = (
                    infsyn.loc[
                        infsyn["parent_infectious_syndrome"].isin(["L3_skin_bacterial"]),
                        "infectious_syndrome",
                    ]
                    .unique()
                    .tolist()
                )
                skin_infect_bacterial += ["L3_skin_bacterial"]

                filtered_df = df[
                    (df["infectious_syndrome"].isin(skin_infect_non_bacterial))
                    & (df[self.int_cause] == 1)
                ]

                if len(filtered_df) == 0:
                    print_log_message("No skin non-bacterial cases to fix")
                initial_deaths = filtered_df.deaths.sum()

                filtered_df[self.int_cause] = filtered_df["all_possible_syndromes"].apply(
                    lambda x: 0 if set(x).isdisjoint(set(skin_infect_bacterial)) else 1
                )

                deaths_removed = filtered_df.loc[filtered_df[self.int_cause] == 0]["deaths"].sum()

                print_log_message(
                    f"Deaths before dropping = {initial_deaths}, unflagged {deaths_removed} non-bacterial skin deaths, {initial_deaths - deaths_removed} remaining (has bacterial)"
                )

                df.update(filtered_df[[self.int_cause]])

        df = df.drop(["_row_id", "l2_syndrome"], axis="columns")

        return df

    def set_part2_flag(self, df):
        p2_cols = [x for x in df.columns if "pII" in x]
        int_cause_chains = [x for x in df.columns if (self.int_cause in x) and ("multiple" in x)]
        p2_chain_dict = dict(list(zip(p2_cols, int_cause_chains)))
        df["pII_" + self.int_cause] = 0
        for p2_col, chain in sorted(p2_chain_dict.items()):
            df.loc[
                (df[chain].isin(self.full_cause_name)) & (df[p2_col] == 1),
                "pII_" + self.int_cause,
            ] = 1
        return df

    def capture_injuries_pattern(self, df, int_cause_cols):
        df["pattern"] = ""
        df["pII_ncodes"] = ""
        df["pII_in_ncodes"] = ""

        p2_cols = [x for x in df.columns if "pII" in x]
        int_cause_chains = [x for x in df.columns if (self.int_cause in x) and ("multiple" in x)]
        p2_chain_dict = dict(list(zip(p2_cols, int_cause_chains)))

        for p2_col, chain in sorted(p2_chain_dict.items()):
            df.loc[df[chain].str.contains("nn", na=False), "pattern"] = (
                df["pattern"] + "_" + df[chain]
            )
            df.loc[(df[chain].str.contains("nn", na=False)) & (df[p2_col] == 1), "pII_ncodes"] = (
                df["pII_ncodes"] + "_" + df[chain]
            )
        df.loc[df["pII_ncodes"] != "", "pII_in_ncodes"] = 1
        return df

    def standard_map_prep(self, df, int_cause_cols):
        bad = df.loc[df["cause_standard_map"] == self.null_code]
        drop_tol = 0.01
        assert len(bad) / len(df["cause_standard_map"]) <= drop_tol
        df = df.loc[df["cause_standard_map"] != self.null_code]
        assert (df.cause_standard_map != self.null_code).all()

        df["standard_map"] = "standard map"
        return df

    def cross_tabulate_shaping(self, df, int_cause_cols):

        print_log_message("Starting cross_tabulate shaping")
        start_deaths = df.deaths.sum()
        assert "ucod" not in df
        assert "_row_id" not in df
        df["_row_id"] = list(range(0, len(df)))

        chain_cols = [
            x
            for x in df.columns
            if ("multiple_cause_" in x) and ("cross_tabulate" in x) and ("pII" not in x)
        ]

        df["ucod_total"] = 5000
        chain_cols += ["ucod_total"]

        id_vars = list(set(df) - set(chain_cols))

        print_log_message("Cross_tabulate: Going from wide to long format")
        df = df.melt(id_vars=id_vars, var_name="position", value_name="int_cause")

        print_log_message("Cross_tabulate: reshape complete, performing ucod tot assertion check")
        ucod_tot = (
            df.groupby(
                by=[
                    "age_group_id",
                    "location_id",
                    "sex_id",
                    "year_id",
                    "cause_cross_tabulate",
                    "int_cause",
                ]
            )["deaths"]
            .sum()
            .reset_index()
        )
        ucod_tot = ucod_tot.loc[ucod_tot["int_cause"] == 5000]

        print_log_message("Cross_tabulate: dropping duplicates")
        df = df.drop_duplicates(subset=["_row_id", "int_cause"], keep="first")
        bad = df.loc[df["cause_cross_tabulate"] == self.null_code]
        drop_tol = 0.01

        print_log_message("Cross_tabulate: assertion checks - ok number of nulls")
        assert len(bad) / len(df["cause_cross_tabulate"]) <= drop_tol
        df = df.loc[df["cause_cross_tabulate"] != self.null_code]
        assert (df.cause_cross_tabulate != self.null_code).all()
        df = df.loc[df.int_cause != self.null_code]

        assert start_deaths == ucod_tot["deaths"].sum()

        df = df.rename(columns={"int_cause": "cross_tabulate"})

        df = df.drop(["_row_id", "position"], axis="columns")
        return df

    def prep_acause_map(self):
        ch = get_current_cause_hierarchy()
        ch = ch[["cause_id", "acause"]]

        mcod_map = ch.set_index("cause_id")["acause"].to_dict()

        return mcod_map

    def crosstab_apply_abcause(self, df):
        gc_df = df.loc[df["cause_id"] == 743]
        nogc_df = df.loc[df["cause_id"] != 743]

        bcause_map = self.prep_int_cause_map()
        bcause_map = MCoDMapper.fill_missing_ICDs(self, bcause_map, "b_cause_id")

        gc_df["cause_cross_tabulate"] = gc_df["cause"].map(bcause_map)
        cause_int_cause_col_dict = MCoDMapper.prep_raw_mapped_cause_dictionary(
            ["cause"], ["cause_cross_tabulate"]
        )

        if self.code_system_id in [1, 6]:
            print_log_message(
                "Trimming ICD codes and remapping underlying cause/primary diagnosis to b_cause -- for crosstabulate"
            )
            gc_df = MCoDMapper.trim_and_remap(self, gc_df, cause_int_cause_col_dict, bcause_map)

        nogc_df["cause_cross_tabulate"] = nogc_df["cause_id"]

        final_df = gc_df.append(nogc_df)
        return final_df

    def remap_to_salmonella_syndromes(self, df, int_cause_col_dict):
        print_log_message("Remapping other salmonella syndromes to parent")
        infsyn = get_infsyn_hierarchy(infsyn_set_version=self.infsyn_set_version)
        other_salmonella_syndromes = infsyn.loc[
            (infsyn["pathogen"].str.contains("salmonella"))
            & ~(infsyn["infsyn_code"].str.startswith("gt", na=False))
            & (infsyn["infectious_syndrome"] != "salmonella_unspecified_unspecific_site"),
            ["infectious_syndrome", "parent_infectious_syndrome"],
        ].drop_duplicates()
        for col in int_cause_col_dict.values():
            df = df.merge(
                other_salmonella_syndromes, how="left", left_on=col, right_on="infectious_syndrome"
            )
            df.loc[df["parent_infectious_syndrome"].notnull(), col] = df[
                "parent_infectious_syndrome"
            ]
            df.drop(columns=["parent_infectious_syndrome", "infectious_syndrome"], inplace=True)
        return df

    def get_computed_dataframe(self, df, map_underlying_cause=True):
        raw_cause_cols = MCoDMapper.get_code_columns(df)
        df = MCoDMapper.fix_icd_codes(df, raw_cause_cols, self.code_system_id)

        if map_underlying_cause:
            print_log_message("Mapping underlying cause/primary diagnosis")
            if self.is_real_code_system:
                cause_map = get_cause_map(
                    code_map_version_id=self.code_map_version_id, **self.cache_options
                )
                code_map = MCoDMapper.prep_cause_map(cause_map)
                code_map = MCoDMapper.fill_missing_ICDs(self, code_map, "code_id")
                df["cause_mapped"] = df["cause"].map(code_map)

                if self.code_system_id in [1, 6]:
                    print_log_message(
                        "Trimming ICD codes and remapping underlying cause/primary diagnosis"
                    )
                    df = MCoDMapper.trim_and_remap(self, df, {"cause": "cause_mapped"}, code_map)

                df = df.rename(columns={"cause_mapped": "code_id"})
                df["code_id"] = df["code_id"].astype(float)
                df = add_code_metadata(
                    df,
                    "cause_id",
                    code_map_version_id=self.code_map_version_id,
                    **self.cache_options,
                )
                report_if_merge_fail(df, "cause_id", "cause")
            else:
                cause_map = pd.read_csv("FILEPATH", encoding="latin1")[["cause", "cause_id"]]
                cause_map = clean_col_strings(cause_map, "cause", deep_clean=False)
                cause_map = cause_map.drop_duplicates()
                df = clean_col_strings(df, "cause", deep_clean=False)
                df = df.merge(cause_map, how="left", on="cause", validate="many_to_one")
                report_if_merge_fail(df, "cause_id", "cause")

        print_log_message("Mapping chain causes")

        if self.int_cause == "cross_tabulate" or self.int_cause == "standard_map":
            df = MCoDMapper.crosstab_apply_abcause(self, df)

            all_cols = MCoDMapper.get_code_columns(df)
            chain_cols_to_map = [x for x in all_cols if "multiple_cause" in x]

            chain_coi_map = self.prep_int_cause_map()
            chain_coi_map = MCoDMapper.fill_missing_ICDs(self, chain_coi_map, "b_cause_id")
            df = MCoDMapper.map_cause_codes(
                self, df, chain_coi_map, self.int_cause, chain_cols_to_map
            )

        else:
            if self.int_cause in self.infectious_syndromes:
                map_type = "infectious_syndrome"
            else:
                map_type = "code_id"
            int_cause_map = self.prep_int_cause_map()
            int_cause_map = MCoDMapper.fill_missing_ICDs(self, int_cause_map, map_type)
            df = MCoDMapper.map_cause_codes(self, df, int_cause_map, self.int_cause)

        int_cause_cols = [x for x in df.columns if self.int_cause in x]
        int_cause_col_dict = MCoDMapper.prep_raw_mapped_cause_dictionary(
            raw_cause_cols, int_cause_cols
        )

        if self.code_system_id in [1, 6]:
            print_log_message("Trimming ICD codes and remapping chain causes")

            if self.int_cause == "cross_tabulate" or self.int_cause == "standard_map":
                if self.int_cause == "cross_tabulate":
                    cause_int_cause_col_dict = MCoDMapper.prep_raw_mapped_cause_dictionary(
                        ["cause"], ["cause_cross_tabulate"]
                    )
                else:
                    cause_int_cause_col_dict = MCoDMapper.prep_raw_mapped_cause_dictionary(
                        ["cause"], ["cause_standard_map"]
                    )

                chain_raw_cause_cols = [x for x in raw_cause_cols if "multiple_cause" in x]
                chain_int_cause_cols = [x for x in int_cause_cols if "multiple_cause" in x]
                chain_int_cause_col_dict = MCoDMapper.prep_raw_mapped_cause_dictionary(
                    chain_raw_cause_cols, chain_int_cause_cols
                )

                df = MCoDMapper.trim_and_remap(self, df, chain_int_cause_col_dict, chain_coi_map)

            else:
                df = MCoDMapper.trim_and_remap(self, df, int_cause_col_dict, int_cause_map)

        if self.int_cause in self.gbd_acauses + self.infectious_syndromes + ["sepsis"]:
            MCoDMapper.assert_all_mapped(df, int_cause_col_dict)
            if self.int_cause in self.infectious_syndromes:
                df = self.remap_to_salmonella_syndromes(df, int_cause_col_dict)

        elif self.int_cause in ["cross_tabulate", "standard_map"]:
            MCoDMapper.assert_all_mapped(df, cause_int_cause_col_dict)
            MCoDMapper.assert_all_mapped(df, chain_int_cause_col_dict)

        if self.int_cause == "sepsis":
            report_if_merge_fail(df, "cause_" + self.int_cause, "cause")

        print_log_message("Remove duplicates")
        cause_cols = [f"cause_{self.int_cause}"] + int_cause_cols
        df = MCoDMapper.drop_duplicated_chain(df, int_cause_cols, self.int_cause)

        print_log_message("Identifying rows with intermediate cause of interest")
        df = self.capture_int_cause(df, int_cause_cols)
        if not self.drop_p2:
            df = self.set_part2_flag(df)

        if self.int_cause != "standard_map":
            drop_cols = ["_".join([x, self.int_cause]) for x in raw_cause_cols if "multiple" in x]

            if self.int_cause == "pathogen" or self.int_cause == "cross_tabulate":
                drop_cols = drop_cols + raw_cause_cols

        else:
            drop_cols = [x for x in raw_cause_cols if x in df.columns]

        check_cols = [x for x in df.columns if "multiple_cause_" in x]
        for col in check_cols:
            if set(df[col]) == set(["0000"]):
                drop_cols += [col]

        drop_cols = list(set(drop_cols))
        drop_cols = [c for c in drop_cols if c in df.columns]

        print_log_message(f"Dropping columns: {drop_cols}")
        df = df.drop(drop_cols, axis=1)

        if self.int_cause == "cross_tabulate":
            int_cols = [x for x in df.columns if x not in ["deaths", "code_id"]]
            for col in int_cols:
                if df[col].dtype != "int64":
                    df[col] = df[col].astype("int64")

        return df
