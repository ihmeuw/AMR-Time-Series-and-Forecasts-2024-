import pandas as pd
import numpy as np
import pyarrow.parquet as pq
import pyarrow as pa
from datetime import datetime
import getpass
user_name = getpass.getuser()
import sys
sys.path.append("ADDRESS")
from claims_repo_clone import claims_config
from function_tracker import tracker
from claims_logger import logger

class process_year:
    def __init__(self, country_config, year=None, sex=None, age=None) -> None:
        self.country = country_config
        self.year = year
        self.sex = sex
        self.age = age
        self.conf = claims_config()
        self.duration_map = pd.read_csv(self.conf.get_path(self.country, "map_file"))
        self.bcause_map = self.duration_map[["code_id",
                                             "min_age","max_age",
                                             "wrong_location","not_disease","sex_restriction",
                                             "b_cause_id","amr_bcause_id"]]
        self.duration_dict = {
            key: value for key, value in np.column_stack((self.duration_map['amr_bcause_id'],
                                                          self.duration_map['duration_day']))}
        self.duration_dict[9999] = 365
        self.code_id_file = "FILEPATH"
        self.function_times = {}
        self.logger = logger(self.country, self.year, self.sex, self.age)
        self.restiction_errors = []
    track_class = tracker()

    @track_class.track_function
    def get_results(self):
        """
        This function loads the long format code_id cleaned parquet file by sex/age
            and provides the outpatient-to-inpatient/outpatient-total ratios by location/year/sex/age/cause.
        Parameters:
        sex (int): 1 for male and 2 for female
        age (int): age group id as defined by AMR, definitions available at age_group_ids.csv
        Returns:
        numpy array: [cause,sex,age,out_to_in,any_outpatient]
        """
        work = self.load_partition(self.code_id_file)
        work = self.reduce_table_to_df(work)
        work = self.format_df_to_int(work)
        work = self.apply_restriction_map(work)
        np_array = self.map_duration_to_numpy(work)
        causes_list = np.unique(np_array[:, 3])
        ratios_list = []
        for c in causes_list:
            cause_arr = np_array[np_array[:, 3] == c]
            if len(cause_arr)==0:
                continue
            n, d, inp, newden = self.calc_outpatientinpatient_event_ratio(cause_arr)
            ratios_list.append([c, self.sex, self.age, n, d, inp, newden])
        return ratios_list

    def load_partition(self, poland_filepath, part=None):
        """
        This function loads the long format code_id cleaned parquet file by sex/age
        Parameters:
        poland_filepath (string): path to parquet file
        sex (int): 1 for male and 2 for female
        age (int): age group id as defined by AMR, definitions available at age_group_ids.csv
        part (int): Default None, if specified will load corresponding 1% of the data in parquet file (for testing)
        Returns:
        pyarrow table : Containing all the data in parquet style format, low memory usage, but very difficult to work with.
        """
        if part != None:
            partitionpath = "FILEPATH"
        else:
            partitionpath = "FILEPATH"
        partition = pq.read_table(partitionpath) 
        table = partition.set_column(16,'code_id', partition['code_id'].cast(pa.int64()))
        return table

    def reduce_table_to_df(self, table, selected=['id', 'date', 'type', 'code_id']):
        """
        This function reduces a pyarrow table to needed columns and converts to pandas dataframe
        Parameters:
        table (pyarrow table): table with all data
        selected (list of strings): Defaults to necessary, different column names can be passed (future flexibility needed?)
        Returns:
        pandas dataframe : Truncated to necessary columns for memory limitation.
        """
        selected_columns = [table[column] for column in selected]
        # Convert the PyArrow arrays to Pandas Series and combine
        selected_series = [selected_column.to_pandas() for selected_column in selected_columns]
        pd_df = pd.concat([pd_s for pd_s in selected_series],axis=1)
        pd_df.columns = selected
        #self.logger.log("info", "Done with reduce_table_to_df")
        return pd_df

    def format_df_to_int(self, pd_df):
        """
        This function reformats the data from objects to int
        Parameters:
        pd_df (pandas DataFrame): required dataframe
        column names could be added for future flexibility?
        Returns:
        pandas dataframe : int datatypes only for increased efficiency and memory usage.
        """
        # Convert from str to timedelta from 1990-01-01 to date
        start_date = datetime(1990, 1, 1)
        dates = pd_df['date']
        pd_df['date'] = pd.to_datetime(dates)
        try:
            pd_df['date'] = (pd_df['date'] - start_date).dt.days
            # Convert the timedelta to integer for easier/faster calculations
            pd_df['date'] = pd_df['date'].astype(int)
        except:
            print(np.unique(pd_df['date']))
            pd_df['date'] = (pd_df['date'] - start_date).dt.days
            pd_df['date'] = pd_df['date'].fillna(0).astype(int)
        pd_df['id'] = pd_df['id'].astype(int)
        pd_df['type'] = pd_df['type'].map({'outpatient':0,'inpatient':1})
        #self.logger.log("info", "Done with format_df_to_int")
        return pd_df.to_numpy().astype(int)

    def apply_restriction(self, c):
        """
        This function applies restiction according to COUNTRY map? and maps code_id to b_cause_id or amr_bcause_id
            *The id for a restricted code_id will map to 9999 and will be correspond to '_miss_code'
        Parameters:
        c (int): amr_bcause_id
        sex (int): 1 for male and 2 for female
        age (int): age group id as defined by AMR, definitions available at FILEPATH
        Returns:
        int: amr_bcause_id
        """
        try:
            row = self.bcause_map.loc[self.bcause_map.code_id == c].iloc[0]
            # Check location restriction
            if row['wrong_location'] == 9:
                return 9999
            # Check disease restriction
            if row['not_disease'] == 9:
                return 9999
            # Check sex restriction
            if row['sex_restriction'] != 3:
                if self.sex != row['sex_restriction']:
                    return 9999
            # Check age restriction
            age_order = [28, 238, 34, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 30, 31, 32, 235]
            index_min_age = age_order.index(row['min_age'])
            index_max_age = age_order.index(row['max_age'])
            index_age = age_order.index(self.age)
            if index_min_age <= index_age <= index_max_age:
                return row['amr_bcause_id']
            else:
                return 9999
        except:
            self.restiction_errors.append(c)
            return 9999

    def apply_restriction_map(self, npdf):
        """
        This function applies restiction according to the dataframe provided.
        Parameters:
        npdf (numpy Array): integer version of data
        Returns:
        numpy Array: Change column 3 amr_bcause_id according to the restriction map.
        """
        code_id_col = npdf[:,3]
        restricted_amrbcombo = np.vectorize(self.apply_restriction)(code_id_col)
        npdf[:,3]= restricted_amrbcombo
        # Fill NaN values with 9999
        nan_indices = np.isnan(npdf[:, 3])
        npdf[nan_indices, 3] = 9999
        #self.logger.log("info", "Done with apply_restriction_map")
        return npdf

    def map_duration_to_numpy(self, np_df):
        """
        This function applies duration according to the map provided.
        Parameters:
        pd_df (pandas DataFrame): integer version of data
        Returns:
        numpy array: numpy array for quick calculations in analysis steps.
        """
        if 822 in self.duration_dict:
            del self.duration_dict[822]
        rows_to_delete = np_df[:, 3] == 822
        np_df = np_df[~rows_to_delete]
        try:
            dur_col = np.vectorize(self.duration_dict.__getitem__)(np_df[:,3])
        except:
            nan_indices = np.isnan(np_df[:,3])
            print("Something going wrong with np.vectorize")
            print(f"Error occured in:\nage: {self.age}\nnan at {nan_indices}")
            for key, value in self.duration_dict.items():
                if np.isnan(value):
                    print(f"Key: {key}, Value: {value}")
            self.logger.log("info", "Something going wrong with np.vectorize")
        try:
            new = np.insert(np_df, 4, dur_col, axis=1)
        except:
            print("Something going wrong with np.insert")
            self.logger.log("info", "Something going wrong with np.insert")
        return new

    def calc_outpatientinpatient_event_ratio(self, cause_arr):
        out_to_in = 0
        any_outpatient = 0
        any_inpatient = 0
        total_events = 0
        for i in np.unique(cause_arr[:,0]):
            filter2 = cause_arr[cause_arr[:,0] == i]
            # If a single instance, no need to check duration
            if len(filter2)==1:
                total_events += 1
                if filter2[0,2]==0:
                    any_outpatient +=1
                if filter2[0,2]==1:
                    any_inpatient += 1
            else:
                # Sort by time
                sorted_indices = np.argsort((filter2[:, 1]))
                sorted_data = filter2[sorted_indices]
                # Calculate time differences between visits and difference in type
                time_diffs = np.diff(sorted_data[:, [1, 2]], axis=0)
                visit_number = 1 # initialize visit number
                visit_type = sorted_data[0,2] # initialize visit type as the first row type value
                
                # Assign visit number and a type{outpatient-only = 0, inpatient-only = 1,outpatient-to-inpatient = 9}
                visit_number_column = [[visit_number, visit_type]] # initialize visit number array
                for time_index in range(0,len(time_diffs)):
                    if time_diffs[time_index,0] > sorted_data[time_index+1,4]: # A new visit outside duration=sorted_data[:,4]
                        visit_number +=1
                        visit_type = sorted_data[time_index+1,2]
                        visit_number_column.append([visit_number, visit_type])
                    else: # Remains in the same visit, within duration window
                        if time_diffs[time_index,1] == 1: # If outpatient to inpatient
                            visit_type=9
                        else: # continue without special flag
                            visit_type= sorted_data[time_index+1,2]
                        visit_number_column.append([visit_number,visit_type])
                # Get all combinations of visits and type status, ie: visit is out-to-in, out-only, in-only or in-to-out
                unique_combinations = np.unique(visit_number_column, axis=0)
                # Count all events with inpatient interaction.
                any_inpatient += len(set(unique_combinations[unique_combinations[:, 1] == 1][:, 0].tolist()))
                # Count number of events with at least one outpatient-to-inpatient instance.
                out_to_in += len(set(unique_combinations[unique_combinations[:, 1] == 9][:, 0].tolist()))
                # Count all events with outpatient interaction.
                any_outpatient += len(set(unique_combinations[(unique_combinations[:, 1] == 0)][:, 0].tolist()))
                total_events += len(set(np.unique(unique_combinations[:, 0])))
        return out_to_in , any_outpatient, any_inpatient, total_events