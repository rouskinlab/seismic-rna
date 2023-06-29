import os

from seismicrna import draw

def beautify_title(title):
    title = title.replace('_', ' ')
    title = title[0].upper() + title[1:]
    return title

def strip_extension(filename):
    return '.'.join(filename.split('.')[:-1])


def write_plot(plot):
    name = strip_extension(plot)
        
    return f"""
.. _{name}:

{beautify_title(name)}
{"-"*len(name)}

{docstring_header(getattr(draw.Study, name))}
                
.. raw:: html
    :file: plots_figs/{plot}
    
.. dropdown:: :fa:`eye,mr-1` **DOCSTRING**: {name}

    .. autofunction:: seismicrna.draw.study.Study.{name}
    

    """
            
def docstring_header(func):
    return func.__doc__.split('\n')[0].strip()


def generate_rst():
    
    with open(gallery_path, 'w') as f:
        f.write("""\nGallery\n=========\n\n\n""")
    
    plots =  os.listdir(os.path.join(os.path.dirname(__file__), 'plots_figs'))
    plots.sort()
    
    for plot in plots:
        
        if not plot.endswith('.html'):
            continue
        
        with open(gallery_path, 'a') as f:
            f.write(write_plot(plot))


def generate_html():

    data = draw.load_dataset()

    study = draw.Study()
    study.df = data
    sample, reference, section, family = study.df.iloc[0][['sample', 'reference', 'section', 'family']]

    path_figs = os.path.join(os.path.dirname(__file__), 'plots_figs')
    
    if not os.path.exists(path_figs):
        os.mkdir(path_figs)
    
    for file in os.listdir(path_figs):
        if file.endswith('.html'):
            os.remove(os.path.join(path_figs, file))

    ################################################################################################
    # Generate HTML plots and save them in the param_docs/source/plots/plots_figs folder
    ################################################################################################
    
    study.mutation_fraction(
        sample = sample,
        reference = reference,
        section='ROI',
        to_html = os.path.join(path_figs, 'mutation_fraction.html')
    )
    
    
    study.mutation_fraction_identity(
        sample = sample,
        reference = reference,
        section='ROI',
        to_html = os.path.join(path_figs, 'mutation_fraction_identity.html')
    )
    
    study.deltaG_vs_sub_rate(
        sample = sample,
        section = 'ROI',
        base_type = ['A','C'],
        to_html = os.path.join(path_figs, 'deltaG_vs_sub_rate.html')
    )
    
    study.experimental_variable_across_samples(
        experimental_variable = 'temperature_k',
        reference = reference,
        section = 'ROI',
        base_type = ['A','C'],
        base_pairing = False,
        to_html = os.path.join(path_figs, 'experimental_variable_across_samples.html'))
    
    # study.auc(
    #     sample = sample,
    #     family = family,
    #     section = 'ROI',
    #     to_html = os.path.join(path_figs, 'auc.html')
    # )
  
    study.mutations_in_barcodes(
        sample = sample,
        to_html = os.path.join(path_figs, 'mutations_in_barcodes.html')
    )
    
    study.num_aligned_reads_per_reference_frequency_distribution(
        sample = sample,
        section = 'full',
        to_html = os.path.join(path_figs, 'num_aligned_reads_per_reference_frequency_distribution.html')
    )
    
    study.mutation_fraction_delta(
        sample1 = study.df['sample'].unique()[0],
        sample2 = study.df['sample'].unique()[1],
        reference1 = reference,
        reference2 = reference,
        section1 = 'ROI',
        section2 = 'ROI',
        to_html = os.path.join(path_figs, 'mutation_fraction_delta.html')
    )

    study.mutations_per_read_per_sample(
        sample = sample,
        section = 'full',
        to_html = os.path.join(path_figs, 'mutations_per_read_per_sample.html')
    )
    
    study.base_coverage(
        sample = sample,
        reference = reference,
        section = 'full',
        to_html = os.path.join(path_figs, 'base_coverage.html')
    )
    
    study.mutation_per_read_per_reference(
        sample = sample,
        reference = reference,
        section = 'full',
        to_html = os.path.join(path_figs, 'mutation_per_read_per_reference.html')
    )
    
    ################################################################################################
        
    
    
def main():
    generate_html()
    generate_rst()

            
if __name__ == '__main__':
    main()
