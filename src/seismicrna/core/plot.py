import sys

CHROME_DEPS_MESSAGE = (
    "Kaleido could not start Chrome for Testing, which SEISMIC-RNA downloads "
    "automatically to export images. On a minimal Linux system (e.g. a "
    "bare-bones Docker image), this usually means a few system libraries "
    "that most systems already have are missing. Install them with:\n"
    "\n"
    "    sudo apt-get install -y libnss3 libnspr4 libatk1.0-0 "
    "libatk-bridge2.0-0 libcups2 libdrm2 libxkbcommon0 libxcomposite1 "
    "libxdamage1 libxfixes3 libxrandr2 libgbm1 libasound2 libpango-1.0-0 "
    "libcairo2 fonts-liberation\n"
    "\n"
    "See the Install page of the SEISMIC-RNA documentation for details."
)


def write_plotly_image(figure, file, **kwargs):
    """Shield sys.argv from Kaleido's logistro dependency, which parses real CLI args at import and can crash (see kaleido==0.2.1 pin history); also ensure Chrome for Testing is downloaded, since Kaleido no longer bundles it."""
    argv = sys.argv
    sys.argv = argv[:1]
    try:
        import kaleido
        from choreographer.errors import BrowserFailedError

        kaleido.get_chrome_sync()
        try:
            figure.write_image(file, **kwargs)
        except BrowserFailedError as error:
            raise RuntimeError(f"{error}\n\n{CHROME_DEPS_MESSAGE}") from error
    finally:
        sys.argv = argv
