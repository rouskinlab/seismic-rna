import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None
import plotly.graph_objects as go

def __find_base_in_sequence(sequence, base_type):
    return [i for i, base in enumerate(sequence) if base in base_type]

def __index_selected(row, base_index, base_type, base_pairing, RNAstructure_use_DMS, RNAstructure_use_temp):
    index = list(range(len(row['sequence'])))
    if base_index is not None:
        if isinstance(base_index, int):
            index = [base_index]
        if isinstance(base_index, list):
            index = base_index
        if isinstance(base_index, str):
            assert (count:=row['sequence'].count(base_index)) == 1, f"{count} sequences {base_index} found in sequence of sample {row['sample']} reference {row['reference']} section {row['section']} cluster {row['cluster']} (sequence: {row['sequence']})"
            temp_idx = row['sequence'].find_files(base_index)
            index = list(range(temp_idx, temp_idx+len(base_index)))

    if base_type is not ['A','C','G','T']:
        index = list(set(index) & set(__find_base_in_sequence(row['sequence'], base_type)))
    
    if base_pairing is not None:
        base_pairing_switch = lambda idx, base, pairing: idx if (base=='.' and not pairing) or (base in ['(',')'] and pairing) else None 
        structure = [base_pairing_switch(i,base,base_pairing) for i, base in enumerate(row['structure'])]
        index = list(set(index) & set([i for i in structure if i is not None]))
    return index


def get_df(df, sample=None, reference=None, section=None, cluster=None, min_cov=0, base_index=None, base_type=['A','C','G','T'], base_pairing=None, 
           RNAstructure_use_DMS=False, RNAstructure_use_temp=False, unique_id = False, index_selected=False, **kwargs)->pd.DataFrame:
    """Get a dataframe with filtered data

    Args:
        df (pd.Dataframe): Dataframe to filter
        sample (list, int, str, optional): Filter rows by sample (list of samples or just a sample). Defaults to None.
        reference (list, int, str, optional): Filter rows by reference (list of references or just a reference). Defaults to None.
        section (list, int, str, optional): Filter rows by section (list of sections or just a section). Defaults to None.
        cluster (list, int, str, optional): Filter rows by cluster (list of clusters or just a cluster). Defaults to None.
        min_cov (int, optional): Filter rows by a minimum threshold for base coverage. Defaults to 0.
        base_index (list, int, str, optional): Filter per-base attributes (sub_rate, sequence, etc) by base index, using 1-indexing. Can be a unique sequence in the row's sequence, a list of indexes or a single index. Gives a  Defaults to None.
        base_type (list, str, optional): Filter per-base attributes (sub_rate, sequence, etc) by base type. Defaults to ['A','C','G','T'].
        base_pairing (bool, optional): Filter per-base attributes (sub_rate, sequence, etc) by expected base pairing. See RNAstructure_use_XXX arguments. Defaults to None.
        RNAstructure_use_DMS (bool, optional): Use DMS for the RNAstructure prediction when filtering by base pairing. Defaults to False.
        RNAstructure_use_temp (bool, optional): Use temperature for the RNAstructure prediction when filtering by base pairing. Defaults to False.
        unique_id (bool, optional): Add a row that will be a unique identifier for each row (combination of sample, reference, section, cluster). Defaults to False.
        index_selected (bool, optional): Add a row that will be a list of indexes selected by the filtering. Defaults to False.
        **kwargs: Additional arguments to pass to filter rows by. Ex: flank='flank_1' will keep only rows with flank=flank_1.

    Returns:
        pd.Dataframe: a filtered dataframe according to the given parameters
    """

    df = df.copy()
    assert df.shape[0] > 0, "Empty dataframe"
    
    if base_index != None:
        if type(base_index) == int:
            base_index = [base_index]
        if hasattr(base_index, '__iter__'):
            base_index = [b-1 for b in base_index] # convert to 0-based index
        else: 
            raise ValueError(f"base_index must be a list, int or None. Got {type(base_index)}")
        assert all([b >= 0 for b in base_index]), f"base_index must be a list of positive integers. Got {base_index}"

    # filter mutation profiles
    if min(df.min_cov) < min_cov:
        df = df.loc[df.min_cov >= min_cov,:]
    for key, value in kwargs.items():
        locals()[key] = value
    mp_attr = ['sample', 'reference', 'section', 'cluster'] + list(kwargs.keys())
    for attr in mp_attr:
        assert attr in df.index, f"Attribute {attr} not found in dataframe"
        if eval(attr) is not None:
            if (isinstance(eval(attr), list) or isinstance(eval(attr), tuple)) or isinstance(eval(attr), np.ndarray):
                df = df[df[attr].isin(eval(attr))]
            else:
                df = df[df[attr] == eval(attr)]
    
    if len(df) == 0:
        return df

        
    # filter base profiles
    if  base_index != None or \
        base_type != ['A','C','G','T'] or \
        base_pairing != None or \
        index_selected:
            
        df['index_selected'] = pd.Series([[]]*df.shape[0], index=df.index)    
        df.loc[:,'index_selected'] = df.apply(lambda row: __index_selected(row, base_index, base_type, base_pairing, RNAstructure_use_DMS, RNAstructure_use_temp), axis=1)
        df = df.loc[df.index_selected.apply(lambda x: len(x) > 0),:]
        bp_attr = ['sequence', 'sub_A','sub_C','sub_G','sub_T', 'sub_N', 'info','del','ins','cov','sub_rate'] + \
            [c for c in df.index.tolist() if (c.startswith('structure') or c.startswith('mod_bases'))]
        for idx, row in df.iterrows():
            for attr in bp_attr:
                # don't filter if the attribute is not an iterable
                if not hasattr(row[attr], '__iter__'):
                    continue
                filtered_cell = [row[attr][i] for i in df.at[idx, 'index_selected']]
                if type(row[attr]) == str:
                    df.at[idx, attr] = ''.join(filtered_cell)
                else:
                    df.at[idx, attr] = np.array(filtered_cell)

    if len(df) == 0:
        return df
    if unique_id:
        try:
            df['unique_id'] = df.apply(lambda row: '_'.join([str(row[attr]) for attr in mp_attr if len(set(df[attr])) > 1]), axis=1)
        except:
            if len(df) > 0:
                "Could not create unique_id column, df: \n {}".format(df)
            else:
                pass
            
    return df
