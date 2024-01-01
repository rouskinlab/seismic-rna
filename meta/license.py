"""

Add a summary of the license to the end of each source file.

"""


import os
from logging import DEBUG, basicConfig, getLogger

basicConfig(level=DEBUG)
logger = getLogger(__name__)


LICENSE_PY = """
########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
"""

LICENSES = {
    ".py": LICENSE_PY,
}


IGNORE_NAME = {".git", ".idea", "__pycache__", "LICENSE"}


def ensure_license(file: str):
    """ Ensure that a source file ends with the license. """
    # Get the license text appropriate for the file type.
    file_name, file_type = os.path.splitext(os.path.basename(file))
    if (license_text := LICENSES.get(file_type)) is None:
        # There is no license text for this file type.
        logger.warning(f"No license text for file {file}")
        return
    # Remove any whitespace surrounding the license text.
    license_text = license_text.strip()
    # Read the file, eliminating trailing whitespace.
    with open(file) as f:
        contents = f.read().rstrip()
    if contents.endswith(license_text):
        # The license already exists.
        logger.debug(f"License already exists in file {file}")
        return
    # Add the license to the end of the file, after one blank line.
    contents = f"{contents}{os.linesep * 2}{license_text}{os.linesep}"
    with open(file, 'w') as f:
        f.write(contents)
    logger.info(f"Added license text to file {file}")


def ensure_licenses(directory: str):
    for name in os.listdir(directory):
        item = os.path.join(directory, name)
        if name in IGNORE_NAME:
            # Ignore the item.
            logger.debug(f"Ignoring {item}")
            continue
        if os.path.isdir(item):
            ensure_licenses(item)
        else:
            ensure_license(item)


def main():
    ensure_licenses(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


if __name__ == "__main__":
    main()

########################################################################
#                                                                      #
# © Copyright 2024, the Rouskin Lab.                                   #
#                                                                      #
# This file is part of SEISMIC-RNA.                                    #
#                                                                      #
# SEISMIC-RNA is free software; you can redistribute it and/or modify  #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation; either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# SEISMIC-RNA is distributed in the hope that it will be useful, but   #
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANT- #
# ABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General     #
# Public License for more details.                                     #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with SEISMIC-RNA; if not, see <https://www.gnu.org/licenses>.  #
#                                                                      #
########################################################################
