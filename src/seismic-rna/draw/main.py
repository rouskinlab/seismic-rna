import os

import pandas as pd
from click import command

from ..core import docdef
from ..core.cli import (opt_out_dir, opt_draw_input,
                        opt_flat, opt_coords, opt_section,
                        opt_mutation_fraction, opt_mutation_fraction_identity,
                        opt_base_coverage, opt_mutations_in_barcodes,
                        opt_mutations_per_read_per_sample,
                        )
from ..core.dump import *
from .study import Study
from logging import getLogger; logger = getLogger(__name__)

params = [
    opt_draw_input,
    opt_out_dir,
    opt_flat,
    opt_coords,
    opt_section,
    opt_mutation_fraction,
    opt_mutation_fraction_identity,
    opt_base_coverage,
    opt_mutations_in_barcodes,
    opt_mutations_per_read_per_sample,
]


@command("draw", params=params)
def cli(**kwargs):
    return run(**kwargs)


@docdef.auto()
def run( 
        inpt: tuple[str, ...],
        *,
        flat: list,
        out_dir: str,
        coords: tuple,
        section: str,
        mutation_fraction: bool,
        mutation_fraction_identity:bool,
        base_coverage: bool,
        mutations_in_barcodes: bool,
        mutation_per_read_per_reference: bool,
        ):
    """Run the draw command.

    """
    if not inpt:
        return
    if type(inpt) == str:
        inpt = (inpt,)
    if type(inpt[0]) == str:
        data = [json.loads(open(i).read()) for i in inpt]  
        logger.info('Loaded {} input files'.format(len(data)))
    else:
        data = inpt
    
    study = Study(
        data = data,
    )
    
    # Check output directory.
    out_dir = os.path.join(out_dir, 'draw')
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
        logger.info('Created output directory: {}'.format(out_dir))
 
    for sample in study.get_samples():
        logger.info('Drawing sample: {}'.format(sample))
            
        path = os.path.join(out_dir, sample)
        if not os.path.isdir(path):
            os.makedirs(path)
            logger.info('Created output directory: {}'.format(path))

        refs, section, base_index = select_coords(coords, sample, study, section)
        
        all_plots = {'mutation_fraction': mutation_fraction, 'mutation_fraction_identity': mutation_fraction_identity, 'base_coverage': base_coverage, 'mutations_in_barcodes': mutations_in_barcodes, 'mutation_per_read_per_reference': mutation_per_read_per_reference}
        logger.info('Drawing {} references - section pairs'.format(len(refs)) +'\n\t'+ '\n\t'.join(['{} - {}'.format(ref, section) for ref, section in zip(refs, section)])\
                    + '\n' + 'plot types: \n\t' + '\n\t'.join([k for k, v in all_plots.items() if v])\
            )

        for ref, section, base_index in zip(refs, section, base_index):
                                
            if flat:
                prefix = os.path.join(path, '__'.join([ref, section,'']))
            else:
                prefix = os.path.join(path, ref, section , '')
                os.makedirs(prefix, exist_ok=True)
                        
            if mutation_fraction:
                study.mutation_fraction(
                    sample=sample,
                    reference=ref,
                    section=section,
                    cluster='pop_avg',
                    #base_index = base_index,
                    to_html = os.path.join(prefix,'mutation_fraction.html') if not flat else prefix + 'mutation_fraction.html',
                )
            
            if mutation_fraction_identity:
                study.mutation_fraction_identity(
                    sample=sample,
                    reference=ref,
                    section=section,
                    cluster='pop_avg',
                    #base_index = base_index,
                    to_html = os.path.join(prefix,'mutation_fraction_identity.html') if not flat else prefix + 'mutation_fraction_identity.html',
                )
            
            if base_coverage:
                study.base_coverage(
                    sample=sample,
                    reference=ref,
                    section=section,
                    cluster='pop_avg',
                    #base_index = base_index,
                    to_html = os.path.join(prefix,'base_coverage.html') if not flat else prefix + 'base_coverage.html',
                )
                
            if mutations_in_barcodes:
                study.mutations_in_barcodes(
                    sample=sample,
                    reference=ref,
                    section=section,
                    cluster='pop_avg',
                    #base_index = base_index,
                    to_html = os.path.join(prefix, 'mutations_in_barcodes.html') if not flat else prefix + 'mutations_in_barcodes.html',
                )
                
            if mutation_per_read_per_reference:
                study.mutation_per_read_per_reference(
                    sample=sample,
                    reference=ref,
                    section=section,
                    cluster='pop_avg',
                    #base_index = base_index,
                    to_html = os.path.join(prefix, 'mutation_per_read_per_reference.html') if not flat else prefix + 'mutation_per_read_per_reference.html',
                )             
                
        logger.info('Done drawing sample: {}'.format(sample))  


def select_coords(coords, sample, study, section):
    
    refs, section_out, base_index = [], [], []

    for t in coords:
        refs.append(t[0])
        section_out.append('full')
        base_index.append(np.arange(t[1], t[2] + 1))
    
    for s in section:
        for ref in study.df[study.df['sample']==sample]['reference'].unique():
            if s in study.df[(study.df['sample']==sample) & (study.df['reference']==ref)]['section'].unique():
                refs.append(ref)
                section_out.append(s)
                base_index.append(None) 

    return refs, section_out, base_index
