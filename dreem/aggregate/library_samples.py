import pandas as pd
import yaml

from dreem.aggregate.resources.get_attributes import read_sample_attributes

def get_samples_info(df_samples, sample, verbose= False):
    if verbose: print(f"Adding samples info for {sample}")

    # Sanity check
    # Keep only the columns that are in sample_attributes.yml
    df_samples = df_samples[df_samples['sample']==sample]

    assert len(df_samples) <= 1, f"{sample} has more than one line in samples.csv"
    assert len(df_samples) == 1, f"{sample} doesn't have a corresponding line in samples.csv"

    exp_env = df_samples['exp_env'][0]
    assert exp_env in ['in_vivo','in_vitro'], f"{exp_env} is not a valid value for exp_env. Should be in_vivo or in_vitro"

    # Load list of mandatory columns and check that they are in samples.csv and not empty for this sample 
    sample_attributes = read_sample_attributes()
    for mand in sample_attributes['mandatory']['all'] + sample_attributes['mandatory'][exp_env]:
        assert mand in list(df_samples.columns), f"{mand} is not in samples.csv"
        assert df_samples[df_samples['sample']==sample][mand].isnull().sum() == 0, f"{mand} is empty in samples.csv for sample {sample}"
    
    df_samples = df_samples.iloc[0]
    
    return df_samples.to_dict()

def get_library_info(df_library, construct, verbose= False):
    if verbose: print(f"Adding library info for {construct}")

    # Sanity check
    df_library = df_library[df_library['construct']==construct]

    df_library.drop(columns = [c for c in df_library.columns if c in ['section_start','section_end','section']], inplace=True)

    for c in df_library.columns:
        assert df_library[c].unique().shape[0] == 1, f"{c} is not unique for construct {construct}"
    
    return df_library.iloc[0].to_dict()
