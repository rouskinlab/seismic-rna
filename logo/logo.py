from PIL import Image

BASE_DIMENSION = 4800
LOGO_FILE = "logo-1440.png"
FAVICON_TEMPLATE = "favicon-{}.ico"
FAVICON_SIZES = [(32, 32)]


def draw_favicon():
    """ Draw the favicon for SEISMIC-RNA. """
    image = Image.open(LOGO_FILE)
    for size in FAVICON_SIZES:
        image.save(FAVICON_TEMPLATE.format("x".join(map(str, size))),
                   format="ICO", sizes=[size])


def draw():
    draw_favicon()


if __name__ == "__main__":
    draw()

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
