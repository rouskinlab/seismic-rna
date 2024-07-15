#!python

"""
Build a Conda package for SEISMIC-RNA.
"""

import os
import re
import shlex
import tomllib
from hashlib import sha256
from urllib.error import URLError
from urllib.request import urlopen

import yaml  # requires pyyaml: pip install pyyaml


def find_script_dir():
    return os.path.dirname(os.path.realpath(os.path.abspath(__file__),
                                            strict=True))


def find_project_dir():
    """ Main directory of SEISMIC-RNA. """
    return os.path.dirname(find_script_dir())


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
    return os.path.join(find_script_dir(), "meta.yaml")


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
            "rnastructure >=6.3",
            "samtools >=1.17"]


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


def write_environment():
    """ Write the environment.yml file for Conda. """
    dependencies: list = list_nonpip_dependencies()
    pip_dependencies = list_pip_dependencies()
    pip_dependencies.append(format_package_version())
    dependencies.append({"pip": list_pip_dependencies()})
    environment = {"name": "seismic",
                   "channels": list_conda_channels(),
                   "dependencies": dependencies}
    with open(find_environment_file(), "w") as f:
        yaml.dump(environment, f)


def format_run_exports_pin():
    return "".join([
        "{{ ", f'pin_subpackage("{find_package_name()}", max_pin="x.x")', " }}"
    ])


def write_metadata():
    """ Write the meta.yaml file for Conda. """
    metadata = {
        "package": {"name": find_package_name(),
                    "version": find_version()},
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
        "channels": list_conda_channels(),
        "test": {"imports": ["seismicrna"]},
        "about": {
            "home": find_github_home(),
            "license": "GPL-3.0",
            "license_family": "GPL",
            "license_file": "LICENSE",
            "license_url": "https://www.gnu.org/licenses/gpl-3.0.html",
            "summary": "SEISMIC-RNA software by the Rouskin Lab"
        },
    }
    with open(find_metadata_file(), "w") as f:
        yaml.dump(metadata, f)


if __name__ == "__main__":
    write_metadata()
