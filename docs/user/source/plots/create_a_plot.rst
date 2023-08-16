
Write your own plots
--------------------

You can write your own plots following these steps.

1. Select your data
********************

Use the ``study.get_df()`` method to get a filtered pandas DataFrame.

.. code::

    data = study.get_df(
            sample = '65degrees_1_S20_L001',        # select one or multiple sample(s)
            reference = ['3042-O-flank_1=hp1-DB',   # select one or multiple reference(s)
                         '3043-CC-flank_1=hp1-DB'],
            section = 'ROI',                        # select one or multiple section(s)
            base_type = ['A','C'],                  # select one or multiple base type(s)
            index_selected = True,                  # add a column with the index of the selected bases (here the bases that are A or C)
        )[['sample','reference','section','sequence','index_selected','sub_rate','deltaG','family','num_aligned','DMS_conc_mM']] # select the columns you want to keep

Output:

 ====================== ======================== ========= ============== ================================= ======================================================================================================================================== ======== ======== ============= ============= 
  sample                 reference                section   sequence        index_selected                   sub_rate                                                                                                                                deltaG   family   num_aligned   DMS_conc_mM  
 ====================== ======================== ========= ============== ================================= ======================================================================================================================================== ======== ======== ============= ============= 
  65degrees_1_S20_L001   3042-O-flank_1=hp1-DB    ROI       AAAACAAAAAAC   [1,2,4,6,7,8,11,12,14,16,17,21]   [0.01309329 0.00612996 0.00246508 0.         0.00243704 0.0174939 0.01204319 0.01748678 0.00246407 0.00367647 0.0020475  0.00493218]     -8.4     hp1      2463          105.0        
  65degrees_1_S20_L001   3043-CC-flank_1=hp1-DB   ROI       AAAACAAAAAAC   [1,2,4,6,7,8,11,12,14,16,17,21]   [0.00958084 0.00646862 0.00240616 0.00168188 0.00502152 0.01409797  0.01364522 0.01192748 0.00506879 0.00239636 0.0021692  0.00621118]   -8.4     hp1      4197          105.0        
 ====================== ======================== ========= ============== ================================= ======================================================================================================================================== ======== ======== ============= ============= 

Take a look at the filtering options in the ``Study.get_df()`` method documentation:

.. dropdown:: :fa:`eye,mr-1` **DOCSTRING**: ``Study.get_df()``

    .. autofunction:: dreem.draw.study.Study.get_df
    
2. Plot your data
******************

Here, we'll use plotly as an example. You can use any plotting library you want.

.. code::

    from plotly import graph_objs as go
    import numpy as np

    def my_plot(data):

        fig = go.Figure()

        # Plot both lines of your DataFrame
        for _, row in data.iterrows():
            fig.add_trace(go.Bar(x=np.arange(0,len(row['sub_rate'])), y=row['sub_rate'], name=row['reference']))

        # Add a title and axis labels
        fig.update_layout(
            title='Comparison of mutation rates between replicates',
            xaxis_title='Base position',
            yaxis_title='Mutation fraction',
        )

        # Show the figure in a Jupyter notebook or in a browser
        fig.show()

        # Write the figure to an HTML file
        fig.write_html('my_figure.html')

        return {'fig': fig, 'data': data}

    my_plot(data)['fig'].show()

You'll get the following output:

.. raw:: html
    :file: my_figure.html


3. Add your plot to DREEM
**************************

This project is community-driven. If you want to add your plot to DREEM, please follow these steps and send us `a pull request <https://github.com/rouskinlab/dreem/pulls>`_.

1. Setup your development environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Fork the DREEM repository on GitHub.
2. Clone your fork locally.
3. Open your favorite terminal and navigate to the cloned repository.
4. Create a new branch for your plot.
5. Install the development dependencies using:

.. code::

    pip install -r requirements.txt

6. In your favorite IDE, open:
    - ``dreem/dreem/draw/study.py`` 
    - ``dreem/dreem/draw/plotter.py`` 
    - ``docs/source/plots/gallery_generator.py``
    - A Jupyter notebook 


2. Example
^^^^^^^^^^^


In this example, we'll add the plot :ref:`mutations_in_barcodes` to DREEM.
You need to add your plot to the following files:

In ``dreem/draw/study.py``:

.. code::

    # In dreem/draw/study.py
    @save_plot
    @doc_inherit(save_plot, style=style_child_takes_over_parent)
    @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
    def base_coverage(self, **kwargs):
        """Plot the base coverage of one or several rows of your dataframe.

        """
        return self.wrap_to_plotter(plotter.base_coverage, locals(), kwargs)

In ``dreem/draw/plotter.py``:

.. code::

    # In dreem/draw/plotter.py
    def base_coverage(data):
        fig = go.Figure()
        [...]
        return {'fig': fig, 'data': data}

In ``docs/source/plots/gallery_generator.py``:

.. code::

    # In docs/source/plots/gallery_generator.py
    def generate_html():
        [...]
        study.base_coverage(
            sample = sample,
            section = 'full',
            to_html = os.path.join(path_figs, 'base_coverage.html'))


3. Add your plot
^^^^^^^^^^^^^^^^^^

1. Add your plot function to ``dreem/draw/plotter.py``.

.. code::

    # In dreem/draw/plotter.py
    def my_plot(data):
        fig = go.Figure()
        [...]
        return {'fig': fig, 'data': data}
    

2. Add your plot to the ``Study`` class in ``dreem/draw/study.py``. 
Use the wrapper: it loads the data for you while making sure that the inputs are valid.

.. code::

    # In dreem/draw/study.py
    class Study:
        [...]
        def my_plot(self, **kwargs):
            """
            My new plot.
            """
            return self.wrap_to_plotter(plotter.my_plot, locals(), kwargs)

3. Add mandatory arguments or default values for optional arguments to your plot function. Document it in the docstring.

.. code::

    # In dreem/draw/study.py
    class Study:
        [...]
        def my_plot(self, sample, reference, section='full', base_type=['A','C'], **kwargs):
            """
            My new plot.

            Args:

                sample (str): Sample name.
                reference (str): Reference name.
                section (str): Section name. Defaults to 'full'.
                base_type (str): Base type. Defaults to ['A','C'].
            """
            return self.wrap_to_plotter(plotter.my_plot, locals(), kwargs)


4. Add the documentation for your plot using the ``custom_dostring.doc_inherit`` decorator. 
When pushing the docs to GitHub Pages, this will add the docstring of the generic plotting function to the docstring of your function.

.. code::

    # In dreem/draw/study.py
    class Study:
        [...]
        # Use this decorator for plots that take one or multiple rows of the DataFrame (use by default).
        @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
        def my_plot(self, sample, reference, section='full', base_type=['A','C'], **kwargs):
        [...]
        # Use this decorator for plots that take a single row of the DataFrame.
        @doc_inherit(default_arguments_single_row, style=style_child_takes_over_parent)
        def my_other_plot(self, sample, reference, section='full', base_type=['A','C'], **kwargs):
        [...]

5. Use the ``@save_plot`` decorator to add the ``to_html`` and ``to_png`` arguments to your plot function.
Add also the documentation for these arguments. Keep the decorators in this order.

.. code::

    # In dreem/draw/study.py
    class Study:
        [...]
        @save_plot
        @doc_inherit(save_plot, style=style_child_takes_over_parent)
        @doc_inherit(default_arguments_multi_rows, style=style_child_takes_over_parent)
        def my_plot(data):
            [...]

6. Test your plot using the testing dataset in a Jupyter notebook.

.. code::

    # In a Jupyter notebook
    from dreem.draw import Study, load_dataset
    study = Study()
    study.df = load_dataset()
    # Plot the first row of the DataFrame
    sample, reference, section, cluster = study.df.iloc[0][['sample', 'reference', 'section', 'cluster']]
    study.my_plot(sample=sample, reference=reference, section=section, cluster=cluster, to_html='my_plot.html')

7. Open ``docs/source/plots/gallery_generator.py``. In ``generate_html()``, generate an HTML file for your plot.

.. code::

    # In docs/source/plots/gallery_generator.py
    def generate_html():
        [...]
        ################################################################################
        # Generate HTML plots and save them in the docs/source/plots/plots_figs folder #
        ################################################################################
        [...]
        study.my_plot(
            sample=sample, 
            reference=reference, 
            section=section, 
            cluster=cluster, 
            to_html='my_plot.html')

8. Run ``docs.source.plots.gallery_generator.py`` to generate the HTML file for your plot. 
Your plot will be ``docs/source/plots/plots_figs/my_plot.html``.
Make sure that it looks good!

9. Make the docs by running the following commands in your terminal:

.. code::

    cd docs
    make html

10. Open ``docs/build/html/index.html``. Your plot should be in the gallery.

11. Commit your changes and push them to GitHub. The docs will be automatically updated on GitHub Pages. Make sure that the docstrings are displayed and that the plot looks good.

12. Send us a pull request to the DREEM repository!
