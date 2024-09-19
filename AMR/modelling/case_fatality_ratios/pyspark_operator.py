import getpass
user_name = getpass.getuser()
import sys
sys.path.append("ADDRESS")
import ast
from claims_repo_clone import claims_config
from claims_repo_clone.poland_ratios_generator import process_year
import pyarrow as pa
from pyspark.sql import SparkSession
import pandas as pd
import numpy as np
from pathlib import Path


class pyspark_operator:
    def __init__(self, country):
        self.conf = claims_config()
        self.country = country
        self.app_name = f"{country}_query_session"
        self.user_dir = Path(f'ADDRESS')
        self.spark = SparkSession.builder.appName(self.app_name).getOrCreate()
        self.inf_syndrome_map = pd.read_csv('FILEPATH')
        self.duration_map = pd.read_csv(self.conf.get_path(self.country, "map_file"))
        self.bcause_map = self.duration_map[["code_id","b_cause_id", "b_cause",
                                             "amr_bcause_id", "amr_b_cause_combo"]]
        self.bcause_dict = self.bcause_map.set_index("amr_bcause_id")["amr_b_cause_combo"].to_dict()
        self.bcause_map = self.duration_map[["code_id","b_cause_id", "b_cause",
                                             "amr_bcause_id", "amr_b_cause_combo"]]
        self.ages = ast.literal_eval(claims_config().get("supplement", "age_group_ids"))
        self.years = ast.literal_eval(claims_config().get(self.country, "years"))
    
    def get_year_data(self, year):
        pq_file = "FILEPATH"
        spdf = self.spark.read.parquet(str(pq_file))
        return spdf
    
    def get_query(self, spark_df, target, target_col):
        conditions = []
        if not isinstance(target, list):
            # Ensure target_col is a legitimate column name
            if target_col not in spark_df.columns:
                raise ValueError(f"The column '{target_col}' is not in this DataFrame.")
            conditions.append(spark_df[target_col] == target)
        else:
            # Check inputs for target(s) and traget column(s)
            if not target or not target_col:
                raise ValueError("Both target and target_col must be lists.")
            for t, col in zip(target, target_col):
                if col not in spark_df.columns:
                    raise ValueError(f"The column '{col}' is not in this DataFrame.")
                # Set the condition based on the pair of target and target_col
                conditions.append(spark_df[col] == t)
        if conditions:
            final_condition = conditions[0]
            if len(conditions)>1:
                for cond in conditions[1:]:
                    final_condition = final_condition & cond
            if spark_df.filter(final_condition).count() > 0:
                filtered_df = spark_df.filter(final_condition)
            else:
                return self.spark.createDataFrame([], schema=spark_df.schema)
        else:
            filtered_df = spark_df
        return filtered_df
    
    def pull_cause_df(self, code_id_list, sex=None, age=None):
        results = []
        years = ast.literal_eval(claims_config().get(self.country, "years"))
        for year in years:
            spark_df = self.get_year_data(year)
            for tv in code_id_list:
                target = [tv]
                target_col = ["code_id"]
                if sex != None:
                    target.extend([sex])
                    target_col.extend(["sex_id"])
                if age != None:
                    target.extend([age])
                    target_col.extend(["age_group_id"])
                rdf = self.get_query(spark_df, target, target_col=target_col)
                #print(f"N rows = {result_df.count()} and N cols = {len(result_df.columns)}")
                pddf = rdf.toPandas()
                results.append(pddf)
        try:
            res_df = pd.concat(results, ignore_index=True)
        except:
            print(f"Empty dataframe for {code_id_list}.")
            res_df = pd.DataFrame()
        return res_df
    
    def get_person_data(self, person_id):
            years = ast.literal_eval(claims_config().get(self.country, "years"))
            alldata = []
            for year in years:
                spark_df = self.get_year_data(year)
                target = [person_id]
                target_col = ['id']
                rdf = self.get_query(spark_df, target, target_col=target_col)
                pddf = rdf.toPandas()
                alldata.append(pddf)
            try:
                res_df = pd.concat(alldata, ignore_index=True)
            except:
                print(f"Empty dataframe for {self.id}.")
                res_df = pd.DataFrame()
            return res_df
    
    def get_individuals_from_year(self, spark_df):
        unique_values = spark_df.select('id').distinct().rdd.flatMap(lambda x: x).collect()
        return unique_values
    
    def get_full_history_ids(self, yearlist):
        full_list = []
        for year in yearlist:
            temp = self.get_year_data(year)
            result = self.get_individuals_from_year(temp)
            print(f"For {year} we have {len(result)} individuals.")
            if year == yearlist[0]:
                full_list = set(result)
            else:
                full_list = full_list.intersection(result)
        return list(full_list)
    
    def get_all_in_out_events(self, yearlist):
        total_inpatient = 0
        total_outpatient = 0
        for year in yearlist:
            df = self.get_year_data(year)
            inpatient_count = df.filter(df.type == 'inpatient').count()
            outpatient_count = df.filter(df.type == 'outpatient').count()
            print(f"For {year} we have {inpatient_count} in-events, and {outpatient_count} out-events.")
            total_inpatient += inpatient_count
            total_outpatient += outpatient_count
        return total_inpatient, total_outpatient

    def convert_to_human_readable(self, df):
        df['amr_bcause_id'] = df['amr_bcause_id'].map(self.bcause_dict)
        df['age_group_id'] = df['age_group_id'].map(self.age_dict)
        df['sex_id'] = df['sex_id'].map({1:"male",2:"female"})
        df = df.rename(columns={'amr_bcause_id': 'b_cause_name', 'age_group_id': 'age_group', 'sex_id': 'sex'})
        return df

    def link_cause_by_years(self, cause, c_col="b_cause_id"):
        if isinstance(cause, (str)):
            code_ids = self.bcause_map.loc[self.bcause_map.amr_b_cause_combo == cause,'code_id'].unique()
        elif isinstance(cause, (int, float)):
            code_ids = self.bcause_map.loc[self.bcause_map.amr_bcause_id == cause,'code_id'].unique()
        else:
            raise ValueError("Either a cause_id int or float, or a amr_b_cause_combo str name must be provided for cause variable.")
        if len(code_ids) == 0:
            print(f"No ICD-10 codes found assoicated with '{cause}'.")
            return None
        ratios_list = []
        for sex in [1,2]:
            for age in self.ages:
                ratio_process = process_year("Poland", sex=sex, age=age)
                tempdf = self.pull_cause_df(code_ids, sex=sex, age=age)
                # Skip any cause/sex/age group combination that is not found in the data.
                if tempdf.empty:
                    continue
                arrow_table = pa.Table.from_pandas(tempdf)
                tempdf = ratio_process.reduce_table_to_df(arrow_table)
                tempdf = ratio_process.format_df_to_int(tempdf)
                tempdf = ratio_process.apply_restriction_map(tempdf)
                array = ratio_process.map_duration_to_numpy(tempdf)
                array[:, 4] = np.where(array[:, 4] >= 365, 3650, array[:, 4])
                # If a given combination is restricted, remove it from the data
                if len(np.unique(array[:,3]))>1:
                    for c_id in np.unique(array[:,3]):
                        temparray = array[array[:, 3] == c_id]
                        out_in, out_all, in_all, total = ratio_process.calc_outpatientinpatient_event_ratio(temparray)
                        ratios_list.append([c_id, cause, sex, age, out_in, out_all, in_all, total])
                else:
                    c_id = np.unique(array[:, 3])[0]
                    out_in, out_all, in_all, total = ratio_process.calc_outpatientinpatient_event_ratio(array)
                    ratios_list.append([c_id, cause, sex, age, out_in, out_all, in_all, total])
        return ratios_list