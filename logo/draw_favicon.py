from PIL import Image

PNG_FILE = "logo-1440.png"
ICO_FILE = "favicon.ico"
SIZE = 32


def draw_favicon():
    """ Draw the favicon for SEISMIC-RNA. """
    image = Image.open(PNG_FILE)
    image.save(ICO_FILE, format="ICO", sizes=[(SIZE, SIZE)])


if __name__ == "__main__":
    draw_favicon()
