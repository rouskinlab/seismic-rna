from subprocess import CalledProcessError

from .shell import run_cmd, WHICH_CMD


def dependency_exists(command: str):
    """ Check whether a dependency exists by running `command --version`
    in the shell, which will succeed only if the dependency exists. """
    try:
        # Try to run 'which command' in the shell.
        run_cmd([WHICH_CMD, command])
    except CalledProcessError:
        # The command failed to run.
        return False
    # The command ran successfully.
    return True


def confirm_dependency(dependency: str, module: str = ""):
    """ Check whether a dependency exists and raise an error if not. """
    if not dependency_exists(dependency):
        by = f"by '{module}' " if module else ""
        raise RuntimeError(f"The dependency '{dependency}' is required {by}but "
                           f"was not found. Please install it (if not yet), "
                           f"place the executable '{dependency}' in your PATH, "
                           f"and ensure the command '{WHICH_CMD} {dependency}' "
                           f"prints the location of the executable.")
