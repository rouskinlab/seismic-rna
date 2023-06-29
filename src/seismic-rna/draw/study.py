from ..draw import manipulator, util, plotter
import pandas as pd
import numpy as np
from ..core.dump import sort_dict, flatten_json
import plotly.graph_objects as go
from custom_inherit import doc_inherit
from ..core.docstring import style_child_takes_over_parent
import os
from .util import save_plot, extract_args
import inspect


class Study(object):
    """A class to store information about a study, i.e a set of samples that are relevant to be studied together.

    Attributes:
        name (str, optional): Short description (<~20 char) of your study. Defaults to None.
        samples (List[str], optional): Names of your study's samples. Defaults to None.
        
    Example:
        >>> study = Study('example',['A1', 'B2', 'B3'])
    """

    attr_list = ['name', 'samples']

    def __init__(self, data=None, min_cov=0, filter_by='sample') -> None:
        """Creates a Study object.

        Args:
            data (dict or list[dict] or pandas.DataFrame, optional): Data to use. Can be a dictionary or list of dictionaries containing DREEM-output jsons, or directly a pandas dataframe. Defaults to None.
            min_cov (int, optional): Minimum number of base coverage for a row to be filtered-in. Defaults to 0.
            filter_by (str, optional): Filter rows by sample or study. When filtered by study, if a row passes the filter, rows with the same 'reference', 'section' and 'cluster' fields for all other samples have a sufficient base coverage. Defaults to 'sample'.            

        Example:
            >>> study = Study(data = {'sample':'mysample',{'reference1': {'section1': {'cluster1': {'sub_N': [100], 'cov': [1000]}}}}},
                              min_cov=1000, 
                              filter_by='sample')
            >>> study = Study(data = (
                                        {'sample':'mysample',{'reference1': {'section1': {'cluster1': {'sub_N': [100], 'cov': [1000]}}}}},
                                        {'sample':'mysample2',{'reference1': {'section1': {'cluster1': {'sub_N': [99], 'cov': [1000]}}}}},
                                        ),
                              min_cov=1000, 
                              filter_by='sample')      
        """
        if data is not None:

            df = pd.DataFrame()

            # If data is a list of json, concatenate them into a single dataframe
            if type(data) is not pd.DataFrame:

                # if data isn't iterable, make it a list
                if not hasattr(data, '__iter__') or isinstance(data, dict):
                    data = [data]

                for sample in data:
                    df = pd.concat([df, pd.DataFrame(flatten_json(sort_dict(sample)))], axis=0)

            # Use the dataframe (loaded or created from json)
            self.set_df(df, min_cov=min_cov, filter_by=filter_by)

        else:
            self.df = None

    def set_df(self, df, min_cov=0, filter_by='sample'):

        self.df = df.reset_index(drop=True)

        self.df = self.df[self.df['min_cov'] >= min_cov]

        if filter_by == 'study':
            self.filter_by_study()

        for attr in ['sample', 'reference']:
            self.df[attr] = self.df[attr].astype(str)

        for attr in ['section', 'cluster']:
            if attr not in self.df.index:
                self.df[attr] = 0

        self.df['deltaG'] = self.df['deltaG'].apply(lambda x: 0.0 if x == 'void' else float(x))

    def filter_by_study(self):
        df = self.df.groupby(['reference', 'section', 'cluster']).filter(
            lambda x: len(self.df['sample'].unique()) == len(x['sample'].unique()))
        self.df = df

    def get_samples(self):
        return self.df['sample'].unique()

    def get_references(self, sample: str):
        return self.df[self.df['sample'] == sample]['reference'].unique()

    def get_sections(self, sample: str, reference: str):
        return self.df[(self.df['sample'] == sample) & (self.df['reference'] == reference)]['section'].unique()

    def get_clusters(self, sample: str, reference: str, section: str):
        return \
        self.df[(self.df['sample'] == sample) & (self.df['reference'] == reference) & (self.df['section'] == section)][
            'cluster'].unique()

    def wrap_to_plotter(self, func, loc, kwargs):

        kwargs = {
            **{k: v for k, v in loc.items() if not k in ['self', 'args', 'kwargs']},
            **kwargs
        }

        """Wrapper for the plot functions."""
        return func(
            manipulator.get_df(self.df,
                               **{k: v for k, v in kwargs.items() if
                                  k in list(self.df.index) + extract_args(manipulator.get_df)}),
            **{k: v for k, v in kwargs.items() if k in extract_args(func)})

    def default_arguments_per_base(self):
        """Default arguments for the plot functions.
        
        Args:
            base_index (list[int], list[str], optional): Filter per-base attributes (sub_rate, sequence, etc) by base index. Can be a unique sequence in the row's sequence, a list of indexes or a single index. Defaults to None.
            base_type (list, str, optional): Filter per-base attributes (sub_rate, sequence, etc) by base type. Defaults to ``['A','C','G','T']``.
            base_pairing (bool, optional): Filter per-base attributes (sub_rate, sequence, etc) by expected base pairing. True will keep only base pairs, False will keep only non-base pairs. Defaults to None.
            **kwargs: Additional arguments to pass to filter rows by. Ex: ``flank='flank_1'`` will keep only rows with ``flank==flank_1``. 

        Returns:
            dict: {``'fig'``: a plotly figure, ``'data'``: a pandas dataframe}
            
        """

    @doc_inherit(default_arguments_per_base, style=style_child_takes_over_parent)
    def default_arguments_single_row(self):
        """Default arguments for the mutiple rows plot functions.
        
        Args:
            sample (str, optional): Selects this sample. Defaults to None.
            reference (str, optional): Selects this reference. Defaults to None.
            section (str, optional): Selects this section. Defaults to ``full``.
            cluster (str, optional): Selects this cluster. Defaults to ``pop_avg``.
            
        """

    @doc_inherit(default_arguments_per_base, style=style_child_takes_over_parent)
    def default_arguments_multi_rows(self):
        """Default arguments for the single row plot functions.
        
        Args:
            sample (list, str, optional): Filter rows by sample (a list of samples or just a sample). Defaults to None.
            reference (list, str, optional): Filter rows by reference (a list of references or just a reference). Defaults to None.
            section (list, str, optional): Filter rows by section (a list of sections or just a section). Defaults to None.
            cluster (list, str, optional): Filter rows by cluster (a list of clusters or just a cluster). Defaults to None.
            
        """

    @doc_inherit(default_arguments_per_base, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def get_df(self, **kwargs):
        """Filter the dataframe by the given arguments."""
        return manipulator.get_df(self.df, **kwargs)

    ############################################################################################################
    # Plot functions                                                                                           #
    ############################################################################################################

    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_single_row, style=style_child_takes_over_parent)
    def mutation_fraction(self, sample, reference, section=None, cluster=None, **kwargs) -> dict:
        """Plot the mutation rates as histograms.

        Args:
            show_ci(bool, optional): Show confidence intervals. Defaults to True.
            
        """
        return self.wrap_to_plotter(
            plotter.mutation_fraction,
            locals(),
            kwargs
        )

    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_single_row, style=style_child_takes_over_parent)
    def mutation_fraction_identity(self, sample, reference, section=None, cluster=None, **kwargs) -> dict:
        """Plot the mutation rates as histograms.

        Args:
            show_ci(bool, optional): Show confidence intervals. Defaults to True.
            
        """
        return self.wrap_to_plotter(
            plotter.mutation_fraction_identity,
            locals(),
            kwargs
        )

    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def deltaG_vs_sub_rate(self, **kwargs) -> dict:
        """Plot the Mutation fraction of each paired-expected base of the ROI for each reference of a sample, w.r.t the deltaG estimation.

        Args:
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].

        """
        return self.wrap_to_plotter(
            plotter.deltaG_vs_sub_rate,
            locals(),
            kwargs
        )

    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def experimental_variable_across_samples(self, experimental_variable, reference, section, **kwargs) -> dict:
        """Plot a given experimental variable vs Mutation fraction across samples for a given reference and section.

        Args:
            experimental_variable (str): Name of the experimental variable to plot.
            models (List[str], optional): Models to fit on the data using scipy.optimize.curve_fit. Under the form ``'lambda x, a, b: a*x+b'`` where ``x`` is the variable. Defaults to [].
            
        """
        index_selected = True
        return self.wrap_to_plotter(
            plotter.experimental_variable_across_samples,
            locals(),
            kwargs
        )

    # @save_plot
    # @doc_inherit(save_plot, style=style_child_takes_over_parent)
    # @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    # def auc(self, **kwargs)->dict:
    #     """Plot the AUC for each mutation profile of the selected data. 

    #     """
    #     unique_id = True
    #     return self.wrap_to_plotter(
    #         plotter.auc,
    #         locals(),
    #         kwargs
    #     )

    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def mutations_in_barcodes(self, sample, section='barcode', **kwargs) -> dict:
        """Plot the number of mutations in the barcode per read of a sample as an histogram.

        """
        return self.wrap_to_plotter(
            plotter.mutations_in_barcodes,
            locals(),
            kwargs
        )

    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def num_aligned_reads_per_reference_frequency_distribution(self, sample, section=None, **kwargs) -> dict:
        """Plot the number of aligned reads per reference as a frequency distribution. x axis is the number of aligned reads per reference, y axis is the count of reference that have this number of aligned reads.

        """

        return self.wrap_to_plotter(
            plotter.num_aligned_reads_per_reference_frequency_distribution,
            locals(),
            kwargs
        )

    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    def mutation_fraction_delta(self, **kwargs) -> dict:
        """Plot the Mutation fraction difference between two mutation profiles.
        
        Args:
            sample1: sample of the first mutation profile.
            sample2: sample of the second mutation profile.
            reference1: reference of the first mutation profile.
            reference2: reference of the second mutation profile.
            section1: section of the first mutation profile.
            section2: section of the second mutation profile.
            cluster1: cluster of the first mutation profile.
            cluster2: cluster of the second mutation profile.
            base_index1: base index of the first mutation profile.
            base_index2: base index of the second mutation profile.
            base_type1: base type of the first mutation profile.
            base_type2: base type of the second mutation profile.
            base_pairing1: base pairing of the first mutation profile.
            base_pairing2: base pairing of the second mutation profile. 
        
        Returns:
            dict: {'fig': a plotly figure, 'data': a pandas dataframe}
        
        """

        df1 = manipulator.get_df(self.df, **{k[:-1]: v for k, v in kwargs.items() if
                                             k.endswith('1') and k[:-1] in list(self.df.index) + list(
                                                 manipulator.get_df.__code__.co_varnames)})
        assert len(df1) > 0, 'No rows found for the first mutation profile.'
        assert len(df1) == 1, 'More than one row found for the first mutation profile.'
        df2 = manipulator.get_df(self.df, **{k[:-1]: v for k, v in kwargs.items() if
                                             k.endswith('2') and k[:-1] in list(self.df.index) + list(
                                                 manipulator.get_df.__code__.co_varnames)})
        assert len(df2) > 0, 'No rows found for the second mutation profile.'
        assert len(df2) == 1, 'Only one row should be selected for the second mutation profile.'
        return plotter.mutation_fraction_delta(pd.concat([df1, df2]).reset_index(drop=True),
                                               **{k: v for k, v in kwargs.items() if
                                                  k in plotter.mutation_fraction_delta.__code__.co_varnames})

    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def mutations_per_read_per_sample(self, sample, section=None, **kwargs) -> dict:
        """Plot the number of mutations per read per sample as an histogram.

        """
        return self.wrap_to_plotter(
            plotter.mutations_per_read_per_sample,
            locals(),
            kwargs
        )

    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_single_row, style=style_child_takes_over_parent)
    def base_coverage(self, sample, reference, section=None, cluster=None, **kwargs):
        """Plot the base coverage of a single row of your dataframe.

        """
        return self.wrap_to_plotter(
            plotter.base_coverage,
            locals(),
            kwargs
        )

    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_single_row, style=style_child_takes_over_parent)
    def mutation_per_read_per_reference(self, sample, reference, section=None, cluster=None, **kwargs) -> dict:
        """Plot the number of mutations per read per reference as an histogram.

        """
        return self.wrap_to_plotter(
            plotter.mutation_per_read_per_reference,
            locals(),
            kwargs
        )
