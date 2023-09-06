from logging import getLogger
from pathlib import Path

from .fqunit import FastqUnit
from ..core.shell import CUTADAPT_CMD, FASTQC_CMD, args_to_cmd, ShellCommand

logger = getLogger(__name__)


def fastqc_cmd(fq_unit: FastqUnit,
               out_prefix: Path, *,
               extract: bool,
               n_procs: int):
    args = [FASTQC_CMD,
            "--threads", n_procs,
            "--extract" if extract else "--noextract",
            "--outdir", out_prefix.parent]
    args.extend(fq_unit.paths.values())
    return args_to_cmd(args)


def run_fastqc(fq_unit: FastqUnit, out_dir: Path, **kwargs):
    step = ShellCommand("running FastQC on", fastqc_cmd)
    return step(fq_unit, out_dir.joinpath(fq_unit.sample), **kwargs)


def cutadapt_cmd(fq_inp: FastqUnit,
                 fq_out: FastqUnit, *,
                 n_procs: int,
                 cut_q1: int,
                 cut_q2: int,
                 cut_g1: str,
                 cut_a1: str,
                 cut_g2: str,
                 cut_a2: str,
                 cut_o: int,
                 cut_e: float,
                 cut_indels: bool,
                 cut_nextseq: bool,
                 cut_discard_trimmed: bool,
                 cut_discard_untrimmed: bool,
                 cut_m: int):
    args = [CUTADAPT_CMD, "--cores", n_procs]
    # Quality trimming
    if cut_nextseq:
        cut_qnext = max(cut_q1, cut_q2)
        if cut_q1 != cut_q2:
            logger.warning("NextSeq trimming takes one quality level, but got "
                           f"two ({cut_q1} and {cut_q2}); using {cut_qnext}")
        args.extend(["--nextseq-trim", cut_qnext])
    else:
        args.extend(["-q", cut_q1])
        if fq_inp.paired:
            args.extend(["-Q", cut_q2])
    # Adapter trimming
    adapters = {"g": cut_g1, "a": cut_a1, "G": cut_g2, "A": cut_a2}
    for arg, adapter in adapters.items():
        if adapter and (fq_inp.paired or arg.islower()):
            for adapt in adapter:
                args.extend([f"-{arg}", adapt])
    args.extend(["-O", cut_o])
    args.extend(["-e", cut_e])
    args.extend(["-m", cut_m])
    if not cut_indels:
        args.append("--no-indels")
    if cut_discard_trimmed:
        args.append("--discard-trimmed")
    if cut_discard_untrimmed:
        args.append("--discard-untrimmed")
    # FASTQ format
    if fq_inp.interleaved:
        args.append("--interleaved")
    # Output files
    output_args = list(zip(("-o", "-p"), fq_out.paths.values(), strict=False))
    for flag, value in output_args:
        args.extend([flag, value])
    # Input files
    args.extend(fq_inp.cutadapt_input_args)
    return args_to_cmd(args)


run_cutadapt = ShellCommand("trimming", cutadapt_cmd)
