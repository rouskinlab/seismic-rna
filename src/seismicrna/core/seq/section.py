import warnings

from .region import Region


class Section(Region):

    def __init__(self, *args, **kwargs):
        warnings.warn(
            "Section is deprecated and will be removed in version 0.25. "
            "Run seismic migrate to convert output files to the new format "
            "and suppress this warning.",
            FutureWarning
        )
        super().__init__(*args, **kwargs)
