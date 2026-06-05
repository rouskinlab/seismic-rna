from pathlib import Path
from typing import Iterable

from ..core.logs import logger


def dry_run(fold_cmds: Iterable[str], parent_dir: Path):
    command_file = parent_dir.joinpath("fold-commands.sh")
    command_text = "".join(f"{cmd}\n" for cmd in fold_cmds)
    with open(command_file, "x") as f:
        f.write(command_text)
    logger.debug("Wrote fold commands to {}:\n{}", command_file, command_text)
    return command_file
