from ..core.seq import parse_fasta
import numpy as np
import pandas as pd
import os


def check_library(library: pd.DataFrame, path_fasta: str, path_save_clean_library: str = None):
    """Sanity check for library.csv
    
    # sections - Check that there is only one value for each reference for each column that's not 'reference',
    'section', 'section_start', 'section_end' - Check that every reference is in the fasta file - Check that there
    are no duplicate rows - If there are multiple sections, copy paste the other attributes so that every row is
    covered - Check that the section_start and section_end are in the correct order - Check that the section_start
    and section_end are within the length of the corresponding sequence in the fasta file. - If there are multiple
    sections, check that the section_start and section_end are unique. If not, remove the duplicates. - If section is
    empty, fill it with the section_start and section_end values separated by an underscore - If a reference has no
    sections, name the section 'full'
    
    # barcodes
    - Check that every barcode is in the fasta file at the given position for the given reference
    - Check that the barcode_start is within the length of the corresponding sequence in the fasta file.
    - Check that the barcode_start is unique for each reference
    - Check that each barcode is unique for each reference
    
    Args:
        library (pandas.DataFrame): content of library.csv
        path_fasta (str): path to the fasta file
        path_save_clean_library (str, optional): path to save the clean library.csv. Defaults to None.
        
    Returns:
        Library (pandas.DataFrame): reframed content of library.csv
    
    """

    a = parse_fasta(path_fasta)
    # fasta = pd.DataFrame({k.decode("utf-8") :v.decode("utf-8")  for k,v in a}, index=[0]).T.reset_index().rename(
    # columns={"index":"reference", 0:"sequence"})
    fasta = pd.DataFrame({k: v.decode("utf-8") for k, v in a}, index=[0]).T.reset_index().rename(
        columns={"index": "reference", 0: "sequence"})

    # Check that the references in the library are in the fasta file
    if not set(library["reference"].unique()).issubset(set(fasta["reference"].unique())):
        raise ValueError("Some references are not in the fasta file")

    # Make sure that barcode_start, section_start and section_end are integers
    cols = []
    cols += ["barcode_start"] if "barcode_start" in library.columns else []
    cols += ["section_start"] if "section_start" in library.columns else []
    cols += ["section_end"] if "section_end" in library.columns else []
    for col in cols:
        library[col] = library[col].apply(lambda x: int(x) if not pd.isnull(x) else x).astype("Int64")

    # Add full section if there is no section or no full
    for ref, sequence in fasta.groupby("reference"):
        if not ref in library["reference"].unique():
            continue
        # assert that no section is named 'full' and the section_start and section_end are not empty or 1-len(sequence)
        if "full" in library.loc[library["reference"] == ref, "section"].unique():
            row = library.loc[(library["reference"] == ref) & (library["section"] == "full")].iloc[0]
            if row["section_start"] == None and row["section_end"] == None:
                row["section_start"] = 1
                row["section_end"] = len(fasta.loc[fasta["reference"] == ref, "sequence"].unique()[0])
                library.loc[(library["reference"] == ref) & (library["section"] == "full"), "section_start"] = row[
                    "section_start"]
                library.loc[(library["reference"] == ref) & (library["section"] == "full"), "section_end"] = row[
                    "section_end"]
            if (row["section_start"] == 1) and (
                    row["section_end"] == len(fasta.loc[fasta["reference"] == ref, "sequence"].unique()[0])).any():
                continue
            raise ValueError(
                f"reference {ref} has a section named 'full' but the section_start and section_end are not empty or "
                f"1-len(sequence)")
        else:
            library = pd.concat([library, pd.DataFrame({"reference": ref, "section": "full", "section_start": 1,
                                                        "section_end": len(
                                                            fasta.loc[fasta["reference"] == ref, "sequence"].unique()[
                                                                0])}, index=[0])], axis=0, ignore_index=True)

    library.reset_index(drop=True, inplace=True)
    # Sections Check that there is only one value for each reference for each column that's not 'reference',
    # 'section', 'section_start', 'section_end'
    if 'section' in library.columns:
        for idx, g in library.groupby("reference"):
            for col in g.columns:
                if col not in ["reference", "section", "section_start", "section_end"]:
                    un = g[col].unique()
                    if np.count_nonzero(~pd.isnull(un)) > 1:
                        raise ValueError(f"{col} has multiple values for {idx}")

        # Check that every reference is in the fasta file
        if not set(library["reference"].unique()).issubset(set(fasta["reference"].unique())):
            raise ValueError("Some references are not in the fasta file")

        # Check that there are no duplicate rows and remove them
        if library.duplicated().any():
            library = library.drop_duplicates()

        # If there are multiple sections, copy paste the other attributes so that every row is covered
        if library["section"].nunique() > 1:
            for idx, g in library.groupby("reference"):
                for col in g.columns:
                    if col not in ["reference", "section", "section_start", "section_end"]:
                        library.loc[library["reference"] == idx, col] = g[col].unique()[0]

    if 'section' in library.columns:
        # Check that the section_start and section_end are in the correct order
        if (library["section_start"] > library["section_end"]).any():
            raise ValueError("section_start is greater than section_end")

        # Check that the section_start is 1 or greater
        if (library["section_start"] < 1).any():
            raise ValueError("section_start is less than 1")

        # Check that section_end is within the length of the corresponding sequence in the fasta file.
        for idx, g in library.groupby("reference"):
            sequence = fasta.loc[fasta["reference"] == idx, "sequence"].unique()[0]
            if (g["section_end"] > len(sequence)).any():
                raise ValueError(
                    f"section_end is greater than the length of the sequence for {idx} and section {max(g['section'].unique())} (sequence length: {len(sequence)})")

        # If there are multiple sections, check that the section_start and section_end are unique. If not, remove the duplicates.
        if library["section"].nunique() > 1:
            if library.duplicated(subset=["reference", "section_start", "section_end"]).any():
                library = library.drop_duplicates(subset=["reference", "section_start", "section_end"])

        # Assert no section is named 'full'
        #        if "full" in library["section"].unique():
        #            raise ValueError("Section cannot be named 'full'")

        # If section, section_start, and section_end are all empty for a certain row, fill section of this row with 'full' and set the section_start and section_end to 0 and the length of the sequence
        for idx, row in library.iterrows():
            if np.count_nonzero(pd.isnull(row[['section', 'section_start', 'section_end']])) == 3:
                library.loc[idx, "section_start"] = 1
                library.loc[idx, "section_end"] = len(
                    fasta.loc[fasta["reference"] == row["reference"], "sequence"].unique()[0])
                library.loc[idx, "section"] = '1-' + str(library.loc[idx, "section_end"])

        # If section is empty but not the section start and end, fill it with the section_start and section_end values separated by an underscore
        if library["section"].isna().any():
            library.loc[library["section"].isna(), "section"] = library.loc[
                                                                    library["section"].isna(), "section_start"].astype(
                str) + "-" + library.loc[library["section"].isna(), "section_end"].astype(str)

    # Barcodes
    if "barcode" in library.columns:
        # Check that the barcode_start is unique for each reference
        if library.groupby("reference")["barcode_start"].nunique().max() > 1:
            raise ValueError("barcode_start is not unique for each reference")

        # Check that the barcode_start is within the length of the corresponding sequence in the fasta file.
        for idx, g in library.groupby("reference"):
            sequence = fasta[fasta["reference"] == idx]["sequence"].unique()[0]
            if (g["barcode_start"] + len(g['barcode']) > len(sequence)).any():
                raise ValueError(
                    f"barcode_start + lenght of barcode is greater than the length of the sequence for {idx}")

        # Check that every barcode is in the fasta file at the given position for the given reference
        for idx, g in library.groupby("reference"):
            barcode = g["barcode"].unique()[0]
            barcode_start = g["barcode_start"].unique()[0] - 1
            if barcode != fasta.loc[fasta["reference"] == idx, "sequence"].unique()[0][
                          barcode_start:barcode_start + len(barcode)]:
                raise ValueError(f"Barcode {barcode} is not in the fasta file at the given position for {idx}")

    # Check taht each attribute is unique for each reference
    for col in library.columns:
        if not col.startswith("section") and not col.startswith("barcode") and col not in ["reference"]:
            if library.groupby("reference")[col].nunique().max() > 1:
                raise ValueError(f"{col} is not unique for each reference")
            for reference in library.reference.unique():
                library.loc[library["reference"] == reference, col] = \
                library.loc[library["reference"] == reference, col].unique()[0]

    if path_save_clean_library is not None and os.path.exists(path_save_clean_library):  #
        library.to_csv(os.path.join(path_save_clean_library, 'clean_library.csv'), index=False)
    return library


def read_sample_attributes():
    return {
        'mandatory': {
            'all': [
                'sample',
                'user',
                'date',
                'exp_env',
                'temperature_k',
                'inc_time_tot_secs',
                'DMS_conc_mM',
            ],
            'in_vitro': [
                'buffer',
            ],
            'in_vivo': [
                'cell_line',
            ],
        },
    }


def check_samples(samples):
    """Check that the samples csv file is correct
    
    - Check that the columns are correct
    - Check that "exp_env" is either "in_vivo" or "in_vitro"
    - Check that the sample names are unique
    - Check that no mandatory columns are missing or empty
    
    Args:
        samples (pandas.DataFrame): The samples csv file
        
    Returns:
        pandas.DataFrame: The samples csv file
    """

    samples_attributes = read_sample_attributes()['mandatory']
    mandatory_columns = samples_attributes['all']
    in_vivo_cols = samples_attributes['in_vivo']
    in_vitro_cols = samples_attributes['in_vitro']

    # Check that the columns are correct
    if not set(mandatory_columns).issubset(set(samples.columns)):
        missing_col = list(set(mandatory_columns) - set(samples.columns))
        raise ValueError(f"Some mandatory columns are missing. Columns are: {missing_col}")

    if 'in_vivo' in samples['exp_env'].unique():
        if not set(in_vivo_cols).issubset(set(samples.columns)):
            missing_col = list(set(in_vivo_cols + mandatory_columns) - set(samples.columns))
            raise ValueError(f"Some mandatory columns are missing. Columns are: {missing_col}")

    if 'in_vitro' in samples['exp_env'].unique():
        if not set(in_vitro_cols).issubset(set(samples.columns)):
            missing_col = list(set(in_vitro_cols + mandatory_columns) - set(samples.columns))
            raise ValueError(f"Some mandatory columns are missing. Columns are: {missing_col}")

    # Check that "exp_env" is either "in_vivo" or "in_vitro"
    if not set(samples['exp_env'].unique()).issubset(set(['in_vivo', 'in_vitro'])):
        raise ValueError("exp_env can only be 'in_vivo' or 'in_vitro'")

    for _, sample in samples.iterrows():
        # Check that the sample names are unique
        if samples[samples['sample'] == sample['sample']].shape[0] > 1:
            raise ValueError(f"Sample name {sample['sample']} is not unique")

        # Check that no mandatory columns are missing or empty
        if sample['exp_env'] == 'in_vivo':
            if sample[in_vivo_cols].isna().any():
                raise ValueError(f"Sample {sample['sample']} has missing or empty mandatory columns")
        elif sample['exp_env'] == 'in_vitro':
            if sample[in_vitro_cols].isna().any():
                raise ValueError(f"Sample {sample['sample']} has missing or empty mandatory columns")
        else:
            raise ValueError(f"Sample {sample['sample']} has an invalid sample_type")

    return samples


def compare_fields(d1, d2, fields):
    assert read_dict(d1, fields) == read_dict(d2,
                                              fields), 'expected and real output are not equal for field {} ({} ≠ {})'.format(
        '/'.join(fields), read_dict(d1, fields), read_dict(d2, fields))


def read_dict(d, fields):
    if fields != []:
        if fields[0] not in d:
            raise ValueError('Key {} not found. Check if it\'s implemented.'.format(fields[0]))
        return read_dict(d[fields[0]], fields[1:])
    return d
