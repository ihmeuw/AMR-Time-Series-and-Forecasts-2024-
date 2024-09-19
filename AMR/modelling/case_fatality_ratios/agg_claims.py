import os
import pandas as pd

def combine_files_to_dataframe(path):
    files = [f for f in os.listdir(path) if '.csv' in f]
    counter = 1
    dfs = []
    for file in files:
        if '_raw_ratios' in file:
            continue
        elif '_done_ratios' in file:
            filepath = os.path.join(path, file)
            df = pd.read_csv(filepath)
            dfs.append(df)
            counter += 1
        else:
            print(file)
    combined_df = pd.concat(dfs, ignore_index=True)
    print(combined_df.b_cause.unique())
    columns_to_convert = combined_df.columns.difference(['b_cause', 'sex', 'age', 'out_to_in_RATIO', 'in_to_all_RATIO'])
    combined_df[columns_to_convert] = combined_df[columns_to_convert].astype(int)
    return combined_df

def process_csv(filepath):
    if '.csv' in filepath:
        df = pd.read_csv(filepath)
        print(f"USED FILE:{filepath}")
    else:
        df = combine_files_to_dataframe(filepath)
        print(f"COMBINED FILES IN:{filepath}")
    try:
        merged_df = df.groupby(['amr_bcause_id', 'b_cause', 'age']).sum(numeric_only=True).reset_index()
    except Exception as e:
        print(f"{df.columns}\nERROR:{e}\n{filepath}")
    merged_df['out_to_in_RATIO'] = merged_df['out_to_in'] / merged_df['ANY_outpatient']
    merged_df['in_to_all_RATIO'] = merged_df['ANY_inpatient'] / merged_df['TOTAL']
    synd = filepath.rsplit('/', 1)[1]
    synd = synd.replace('_allyear_ratios.csv', '')
    synd_df = df.groupby(['age']).sum().reset_index()
    synd_df['b_cause'] = synd
    synd_df['amr_bcause_id'] = 9000
    # Calculate the ratio columns
    synd_df['out_to_in_RATIO'] = synd_df['out_to_in'] / synd_df['ANY_outpatient']
    synd_df['in_to_all_RATIO'] = synd_df['ANY_inpatient'] / synd_df['TOTAL']
    result_df = merged_df.append(synd_df, ignore_index=True)
    return result_df

folder = "ADDRESS"

donelist = []
all_df = []
file_folder_list = os.listdir(folder)
for filename in file_folder_list:
    if 'blood_stream' in filename:
        continue
    if 'raw_ratios' in filename:
        continue
    if filename.endswith("FILEPATH"):
        sname = filename.split('_', 1)[1]
        sname = sname.replace('FILEPATH', '')
        results = process_csv(os.path.join(folder, filename))
        results['synd'] = sname
        all_df.append(results)
        donelist.append(sname)
        print(f"{sname} summary saved successfully.")
    elif os.path.isdir(os.path.join(folder, filename)):
        sname = filename.split('_', 1)[1]
        if sname in donelist:
            continue
        results = process_csv(os.path.join(folder, filename))
        results['synd'] = sname
        print(f"{sname} summary saved successfully.")
        all_df.append(results)
    else:
        print("Unknown file -> ",filename)

result = pd.concat(all_df, ignore_index=True)
synd_only = result.loc[result.amr_bcause_id == 9000]
save_path = os.path.join(folder, 'FILEPATH')
synd_only.to_csv(save_path, index=False)
print(f"DONE\nSAVE TO:\n{save_path}")