from setuptools import setup, Extension

module = Extension(
    "relate",  # Module name
    sources=["relate.c"],  # C source file(s)
)

setup(
    name="relate",
    version="1.0",
    description="Relate C extension",
    ext_modules=[module],
)

