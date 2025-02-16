from shutil import which
from os import environ


class DependencyError(RuntimeError):
    """ A required dependency is missing. """


def dependency_exists(dependency: str) -> bool:
    """ Check whether a dependency exists. """
    return which(dependency) is not None


def require_dependency(dependency: str, module: str = ""):
    """ If a dependency does not exist, return an error message. """
    if not dependency_exists(dependency):
        by = f"by '{module}' " if module else ""
        message = (f"{repr(dependency)} is required {by}but was not found. "
                   f"Please install it (if not yet) and place the executable "
                   f"for {repr(dependency)} in your PATH.")
        raise DependencyError(message)

def require_env_var(dependency: str, module: str = ""):
    """ If a required env var is not set, return an error message. """
    if not environ.get(dependency):
        by = f"by '{module}' " if module else ""
        message = (f"The environmental variable {repr(dependency)} is "
                   f"required {by}but was not "
                   "found. Please set it to the "
                   "correct value before continuing."
                   if dependency != "RNARTISTCORE" else
                   f"To use the {module} module, please install RNArtistCore "
                   "manually. To do so, download the latest RNArtistCore .jar "
                   f"file from https://github.com/fjossinet/RNArtistCore/releases "
                   f"and set the environmental variable RNARTISTCORE to the .jar "
                   "path by running: "
                   'export RNARTISTCORE="/full/path/to/your/'
                   'rnartistcore-X.X.X-SNAPSHOT-jar-with-dependencies.jar" '
                   "It is recommended that you add this command to your "
                   ".bashrc or .zshrc file so you don't see this error again. "
                   "If you continue to experience issues, ensure you have "
                   "Java Development Kit (JDK) >= 11 installed. You can check by "
                   "running java --version in your terminal. If problems "
                   "persist, please open an issue at "
                   "https://github.com/rouskinlab/seismic-rna")
        raise DependencyError(message)
