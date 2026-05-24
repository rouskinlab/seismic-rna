from inspect import getmembers

from click import Argument, Option

from . import cli

# Get every parameter defined for the command line interface.
cli_args = dict(getmembers(cli, lambda member: isinstance(member, Argument)))
cli_opts = dict(getmembers(cli, lambda member: isinstance(member, Option)))

# Get the default value for every parameter.
cli_defaults = {param.name: param.default
                for param in (cli_args | cli_opts).values()
                if param.default is not None}

defaults_to_none = dict(fastp_adapter_fasta=None,
                        regions_file=None,
                        min_mut_gap=None,
                        min_mut_gap_weights=None,
                        injected_mut_probs=None,
                        mask_a=None,
                        mask_c=None,
                        mask_g=None,
                        mask_u=None,
                        join_clusts=None,
                        fold_regions_file=None,
                        fold_constraint=None,
                        fold_commands=None,
                        eddy_prior_paired_file=None,
                        eddy_prior_unpaired_file=None,
                        samples_meta=None,
                        refs_meta=None,
                        read_pos=None,
                        collate_out_dir=None,
						seed=None)
