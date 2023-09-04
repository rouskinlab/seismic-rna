"""

Alignment Initialization Module

========================================================================

Expose cli, params, and run at the level of the `align` sub-package, so
that they can be imported in any of the following manners:

>>> import seismicrna
>>> seismicrna.align.run

or

>>> from seismicrna import align
>>> align.run

or

>>> from seismicrna.align import run
>>> run

------------------------------------------------------------------------

"""

from . import main, fq2bam, fqops
from .main import cli, params, run
