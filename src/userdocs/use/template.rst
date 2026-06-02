:orphan:

..
    This file is the reference template for SEISMIC-RNA command-reference
    pages under ``use/workflow/`` and ``use/utility/``.
    It is marked ``:orphan:`` so it does not appear in any toctree.
    When adding or rewriting a command page, mirror this section order.

********************************************************************************
seismic <command>
********************************************************************************

Purpose
================================================================================

One paragraph: what this command does and where it sits in the pipeline.
Mention any preceding/following step it depends on.


Inputs
================================================================================

List every required input file type, where it comes from (which earlier
command produced it, or that it is user-provided), and cross-link the file
format spec under :doc:`/formats/index`.


Outputs
================================================================================

List every output file the command writes, including its on-disk path
template (e.g. ``{out}/{sample}/<cmd>/{ref}/...``) and report file.
Cross-link the relevant format pages under :doc:`/formats/report/index`.


Quick example
================================================================================

A minimal, copy-paste-able invocation that demonstrates the command's
typical use::

    seismic <command> [common opts] <inputs>


Options
================================================================================

Every option, grouped (Input/Output, Filtering, Performance, etc.).
For each option, give the meaning, the default, when to change it, and any
interactions or conflicts with other options.
Long option lists may be split into subsections.


Caveats
================================================================================

Assumptions, ordering requirements, units, and edge-case behavior that a
new user might miss.
Note any renames between releases (cross-link the project ``CHANGELOG.md``).


Performance tips
================================================================================

Parallelism (``--num-cpus``), batching (``--batch-size``), memory
considerations, and other tuning advice.


Common errors
================================================================================

A table or list of exception messages produced by this command, the
underlying cause, and how to fix it.
Mine from the module's ``raise`` and ``class .*Error`` definitions.


Common unexpected results
================================================================================

Outcomes that are not errors but may surprise the user: empty output,
all reads filtered out, NaNs in tables, single-cluster collapse, etc.
Each entry: what the user sees, why it happens, what to do.


See also
================================================================================

Cross-links to the preceding and following pipeline commands, the relevant
file format pages, and the algorithm pages.
