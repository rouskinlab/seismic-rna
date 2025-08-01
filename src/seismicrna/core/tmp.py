import errno

from functools import wraps
from inspect import Parameter, Signature
from pathlib import Path
from shutil import move
from typing import Callable

from .logs import logger
from .path import (mkdir_if_needed,
                   rmdir_if_needed,
                   randdir,
                   sanitize,
                   transpath)

PENDING = "release"
WORKING = "working"


def release_to_out(out_dir: Path,
                   release_dir: Path,
                   initial_path: Path):
    """ Move temporary path(s) to the output directory. """
    logger.routine(
        f"Began releasing {initial_path} from {release_dir} to {out_dir}"
    )
    # Determine the path in the output directory.
    out_path = transpath(out_dir, release_dir, initial_path)
    if initial_path.exists():
        # Ensure the parent directory of the new path exists.
        mkdir_if_needed(out_path.parent)
        # If the output path already exists, then first rename it.
        delete_path = randdir(out_path.parent, f"{out_path.name}-")
        try:
            out_path.rename(delete_path)
        except FileNotFoundError:
            # The output path does not yet exist.
            deleted = False
            logger.detail(f"Output path {out_path} does not yet exist")
        else:
            deleted = True
            logger.action(
                f"Moved output path {out_path} to {delete_path} (to be deleted)"
            )
        try:
            # Move the initial path to the output location.
            try:
                initial_path.rename(out_path)
            except OSError as e:
                if e.errno == errno.EXDEV:
                    logger.warning(
                        "Non-atomic move: temporary and output are on "
                        f"different filesystems; moving {initial_path} to "
                        f"{out_path} via copy-delete; data loss may occur "
                        "if interrupted"
                    )
                    move(initial_path, out_path)
                else:
                    raise
            logger.action(
                f"Moved initial path {initial_path} to output path {out_path}"
            )
        except Exception:
            if deleted:
                # If an error occurred, then restore the original output
                # path before raising the exception.
                delete_path.rename(out_path)
                logger.action(f"Moved {delete_path} (to be deleted) "
                              f"back to output path {out_path}")
            else:
                # No original files were moved to the delete directory,
                # which is therefore still empty. Delete it.
                rmdir_if_needed(delete_path)
            raise
        # Once the initial path has been moved to its destination, the
        # original directory can be deleted safely.
        rmdir_if_needed(delete_path,
                        rmtree=True,
                        raise_on_rmtree_error=False)
    else:
        logger.detail(f"Skipped releasing {initial_path} (does not exist)")
    if not out_path.exists():
        raise FileNotFoundError(out_path)
    logger.routine(f"Ended releasing {initial_path} to {out_path}")
    return out_path


def get_release_working_dirs(tmp_dir: Path):
    release_dir = tmp_dir.joinpath(PENDING)
    working_dir = tmp_dir.joinpath(WORKING)
    release_dir.mkdir(parents=False, exist_ok=True)
    working_dir.mkdir(parents=False, exist_ok=True)
    return release_dir, working_dir


def with_tmp_dir(pass_keep_tmp: bool):
    """ Make a temporary directory, and delete it after returning. """

    def decorator(func: Callable):

        @wraps(func)
        def wrapper(*args,
                    tmp_pfx: str | Path,
                    keep_tmp: bool,
                    **kwargs):
            tmp_dir = None
            try:
                tmp_pfx = sanitize(tmp_pfx)
                tmp_dir = randdir(tmp_pfx.parent, prefix=tmp_pfx.name)
                if pass_keep_tmp:
                    kwargs = dict(keep_tmp=keep_tmp, **kwargs)
                return func(*args, tmp_dir=tmp_dir, **kwargs)
            finally:
                if tmp_dir is not None and not keep_tmp:
                    rmdir_if_needed(tmp_dir,
                                    rmtree=True,
                                    raise_on_rmtree_error=False)

        # Add tmp_pfx and keep_tmp to the signature of the wrapper, and
        # remove tmp_dir (functools.wraps does not do so automatically).
        params = dict(Signature.from_callable(func).parameters)
        params.pop("tmp_dir")
        for name in ["tmp_pfx", "keep_tmp"]:
            if name not in params:
                params[name] = Parameter(name, Parameter.KEYWORD_ONLY)
        # Ensure the variadic keyword parameter (if it exists) is last.
        kwargs_name = None
        for name, param in params.items():
            if param.kind == param.VAR_KEYWORD:
                assert kwargs_name is None
                kwargs_name = name
        if kwargs_name is not None:
            params[kwargs_name] = params.pop(kwargs_name)
        wrapper.__signature__ = Signature(parameters=list(params.values()))
        return wrapper

    return decorator
