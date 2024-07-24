#!python

"""
Build a Conda package for SEISMIC-RNA.
"""

import os
import re
import tomllib
from hashlib import sha256
from shutil import copy2
from urllib.error import URLError
from urllib.request import urlopen


def mkdir_if_needed(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        pass


def find_project_dir():
    """ Main directory of SEISMIC-RNA. """
    return os.path.dirname(os.path.realpath(os.path.abspath(__file__),
                                            strict=True))


def find_conda_dir():
    """ Directory for Conda recipe files. """
    return os.path.join(find_project_dir(), "conda")


def find_git_dir():
    """ Directory for all GitHub repositories. """
    return os.path.dirname(find_project_dir())


def find_bioconda_recipes_dir():
    """ Directory for the bioconda-recipes GitHub repository. """
    return os.path.join(find_git_dir(), "bioconda-recipes")


def find_bioconda_recipe_dir():
    """ Directory for this project's Bioconda recipe. """
    return os.path.join(find_bioconda_recipes_dir(),
                        "recipes",
                        find_package_name())


def find_package_name():
    return os.path.basename(find_project_dir())


def find_pyproject_file():
    """ Pyproject file. """
    return os.path.join(find_project_dir(), "pyproject.toml")


def find_environment_file():
    """ Environment specification file. """
    return os.path.join(find_project_dir(), "environment.yml")


def find_metadata_file():
    """ Build metadata file. """
    return os.path.join(find_conda_dir(), "meta.yaml")


def find_build_script():
    """ Build metadata file. """
    return os.path.join(find_conda_dir(), "build.sh")


def find_version_file():
    """ File containing version of SEISMIC-RNA in this directory. """
    return os.path.join(find_project_dir(),
                        "src",
                        "seismicrna",
                        "core",
                        "version.py")


def find_version():
    """ Determine the version of SEISMIC-RNA in this directory. """
    # Parse the source file instead of importing seismicrna because
    # seismicrna is not guaranteed to be installed, and if it is then
    # the installed version is not necessarily the one being built here.
    pattern = re.compile('__version__ = "([0-9a-z.]+)"')
    version = ""
    version_file = find_version_file()
    with open(version_file) as f:
        for line in f:
            match = pattern.match(line)
            if match:
                if version:
                    raise ValueError(
                        f"Version is defined more than once in {version_file}"
                    )
                version = match.groups()[0]
    if not version:
        raise ValueError(f"Version is not defined in {version_file}")
    return version


def format_package_version():
    """ Format the name and version of this package. """
    return f"{find_package_name()}={find_version()}"


def list_pip_dependencies():
    """ List the pip dependencies in the pyproject.toml file. """
    with open(find_pyproject_file(), "rb") as f:
        return tomllib.load(f)["project"]["dependencies"]


def list_nonpip_dependencies():
    """ List the dependencies not in the pyproject.toml file. """
    return ["python >=3.10",
            "bowtie2 >=2.5.1",
            "fastqc >=0.12.1",
            "rnastructure >=6.2",
            "samtools >=1.17",
            "brotli-python >=1.0"]


def list_all_dependencies():
    return list_nonpip_dependencies() + list_pip_dependencies()


def list_conda_channels():
    return ["bioconda", "conda-forge"]


def find_github_home():
    return f"https://github.com/rouskinlab/{find_package_name()}"


def find_github_file():
    return f"{find_github_home()}/archive/refs/tags/v{find_version()}.tar.gz"


def calc_github_file_sha256():
    url = find_github_file()
    response = urlopen(url)
    if response.status != 200:
        raise URLError(f"{url} returned status {response.status}")
    return sha256(response.read()).hexdigest()


def _is_listlike(item):
    return isinstance(item, (list, tuple, set))


def _is_dictlike(item):
    return isinstance(item, dict)


def _is_iterlike(item):
    return _is_listlike(item) or _is_dictlike(item)


def _make_prefix(indent: int, is_list: bool):
    prefix = "  " * indent
    if is_list:
        prefix += "- "
    return prefix


def _generate_yaml_lines(data: dict | list | set | tuple, indent: int):
    if _is_dictlike(data):
        prefix = _make_prefix(indent, False)
        for key, value in data.items():
            if _is_iterlike(value):
                yield f"{prefix}{key}:"
                yield from _generate_yaml_lines(value, indent + 1)
            else:
                yield f"{prefix}{key}: {value}"
    elif _is_listlike(data):
        prefix = _make_prefix(indent, True)
        for value in data:
            yield f"{prefix}{value}"
    else:
        raise TypeError(data)


def format_yaml_text(data: dict):
    """ Format a dictionary into YAML text. """
    if not isinstance(data, dict):
        raise TypeError(data)
    return "\n".join(["---", ""] + list(_generate_yaml_lines(data, 0)))


def write_environment():
    """ Write the environment.yml file for Conda. """
    environment = {"name": "seismic",
                   "channels": list_conda_channels(),
                   "dependencies": list_nonpip_dependencies()}
    yaml_text = format_yaml_text(environment)
    with open(find_environment_file(), "w") as f:
        f.write(yaml_text)


def format_run_exports_pin():
    return "".join([
        "{{ ", f'pin_subpackage("{find_package_name()}", max_pin="x.x")', " }}"
    ])


def write_metadata():
    """ Write the meta.yaml file for Conda. """
    metadata = {
        "package": {"name": find_package_name(),
                    "version": find_version()},
        "about": {
            "home": find_github_home(),
            "license": "GPL-3.0-only",
            "license_family": "GPL3",
            "license_file": "LICENSE",
            "license_url": "https://www.gnu.org/licenses/gpl-3.0.html",
            "summary": "SEISMIC-RNA software by the Rouskin Lab"
        },
        "source": {"url": find_github_file(),
                   "sha256": calc_github_file_sha256()},
        "build": {
            "noarch": "python",
            "number": 0,
            "run_exports": [format_run_exports_pin()]
        },
        "requirements": {"build": ["python >=3.10",
                                   "hatch >=1.12"],
                         "run": list_all_dependencies()},
        "test": {"imports": ["seismicrna"]},
    }
    yaml_text = format_yaml_text(metadata)
    mkdir_if_needed(find_conda_dir())
    with open(find_metadata_file(), "w") as f:
        f.write(yaml_text)


def write_build_script():
    """ Write the build.sh file for Conda. """
    build_script = "\n".join([
        "#!/bin/bash",
        "",
        "# DO NOT RUN THIS SCRIPT YOURSELF!",
        "# It should only be run by conda build.",
        "",
        "set -euxo pipefail",
        "",
        "$PYTHON -m pip install --no-dependencies $PWD"
    ])
    mkdir_if_needed(find_conda_dir())
    with open(find_build_script(), "w") as f:
        f.write(build_script)


def copy_recipe_to_bioconda():
    """ Copy the recipe files to the Bioconda recipe directory. """
    src = find_conda_dir()
    dst = find_bioconda_recipe_dir()
    mkdir_if_needed(dst)
    for file in os.listdir(src):
        copy2(os.path.join(src, file), dst)


if __name__ == "__main__":
    write_environment()
    write_metadata()
    write_build_script()
    copy_recipe_to_bioconda()
