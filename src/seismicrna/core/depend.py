from shutil import which


def dependency_exists(dependency: str) -> bool:
    """ Check whether a dependency exists. """
    return which(dependency) is not None


def require_dependency(dependency: str, module: str = ""):
    """ Check whether a dependency exists and raise an error if not. """
    if not dependency_exists(dependency):
        by = f"by '{module}' " if module else ""
        raise RuntimeError(f"The dependency '{dependency}' is required {by}but "
                           f"was not found. Please install it (if not yet) and "
                           f"place the executable '{dependency}' in your PATH.")
