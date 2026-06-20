#!/usr/bin/env python3
"""
Move module-level `import numpy as np` and `import pandas as pd` into the
function/method bodies that use them. Also adds `from __future__ import
annotations` so that type-hint references to numpy/pandas types (e.g.
np.ndarray, pd.Series) remain valid without a module-level import.

Usage: python defer_np_pd.py [src_dir]
"""

import ast
import re
import sys
from pathlib import Path

NP_RE = re.compile(r"\bnp\b")
PD_RE = re.compile(r"\bpd\b")


def is_docstring(node: ast.stmt) -> bool:
    return (
        isinstance(node, ast.Expr)
        and isinstance(node.value, ast.Constant)
        and isinstance(node.value.value, str)
    )


def collect_functions(nodes: list[ast.stmt], out: list[ast.FunctionDef]) -> None:
    """Collect top-level and class-body FunctionDefs; do not recurse into
    nested function bodies (inner functions close over the outer import)."""
    for node in nodes:
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            out.append(node)
        elif isinstance(node, ast.ClassDef):
            collect_functions(node.body, out)


def body_text(lines: list[str], func: ast.FunctionDef) -> str:
    start = func.body[0].lineno - 1  # 0-indexed
    end = func.body[-1].end_lineno  # exclusive upper bound (0-indexed)
    return "".join(lines[start:end])


def has_module_level_usage(tree: ast.Module, check_np: bool, check_pd: bool) -> bool:
    """Return True if np/pd names appear at module level outside imports,
    functions, and classes (e.g. a module-level constant using numpy)."""
    for node in tree.body:
        if isinstance(
            node,
            (
                ast.Import,
                ast.ImportFrom,
                ast.FunctionDef,
                ast.AsyncFunctionDef,
                ast.ClassDef,
            ),
        ):
            continue
        for child in ast.walk(node):
            if isinstance(child, ast.Name):
                if check_np and child.id == "np":
                    return True
                if check_pd and child.id == "pd":
                    return True
    return False


def transform(path: Path) -> str | None:
    """
    Transform one file in-place.
    Returns None if unchanged, "CHANGED" if modified, or a "SKIP …" message.
    """
    source = path.read_text(encoding="utf-8")
    lines = source.splitlines(keepends=True)

    try:
        tree = ast.parse(source)
    except SyntaxError as e:
        return f"SKIP (syntax error: {e})"

    # ── locate module-level numpy / pandas imports ────────────────────────
    check_np = False
    check_pd = False
    import_line_indices: set[int] = set()  # 0-indexed

    for node in tree.body:
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == "numpy":
                    check_np = True
                    import_line_indices.add(node.lineno - 1)
                elif alias.name == "pandas":
                    check_pd = True
                    import_line_indices.add(node.lineno - 1)

    if not check_np and not check_pd:
        return None

    # ── guard: skip files with module-level np/pd constants ───────────────
    if has_module_level_usage(tree, check_np, check_pd):
        return "SKIP (module-level numpy/pandas usage — needs manual refactor)"

    # ── find __future__ annotations status and module docstring ───────────
    has_future_ann = False
    last_future_end = -1  # 0-indexed end line of last from __future__ …
    module_doc_end = -1  # 0-indexed end line of module docstring

    if tree.body and is_docstring(tree.body[0]) and tree.body[0].col_offset == 0:
        module_doc_end = tree.body[0].end_lineno - 1

    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module == "__future__":
            last_future_end = node.end_lineno - 1
            for alias in node.names:
                if alias.name == "annotations":
                    has_future_ann = True

    # ── collect function definitions (not nested) ─────────────────────────
    functions: list[ast.FunctionDef] = []
    collect_functions(tree.body, functions)

    # ── build edit maps ───────────────────────────────────────────────────
    # inserts_before[i] = list of lines to emit before 0-indexed line i
    inserts_before: dict[int, list[str]] = {}
    lines_to_delete: set[int] = set(import_line_indices)

    def prepend_insert(idx: int, content: list[str]) -> None:
        inserts_before.setdefault(idx, [])
        inserts_before[idx] = content + inserts_before[idx]

    # add from __future__ import annotations
    if not has_future_ann:
        if last_future_end >= 0:
            future_idx = last_future_end + 1
        elif module_doc_end >= 0:
            future_idx = module_doc_end + 1
        else:
            future_idx = 0
        prepend_insert(future_idx, ["from __future__ import annotations\n"])

    # add per-function imports
    for func in functions:
        if not func.body:
            continue

        txt = body_text(lines, func)
        needs_np = check_np and bool(NP_RE.search(txt))
        needs_pd = check_pd and bool(PD_RE.search(txt))
        if not needs_np and not needs_pd:
            continue

        # skip if already imported inside this function
        if needs_np and "import numpy" in txt:
            needs_np = False
        if needs_pd and "import pandas" in txt:
            needs_pd = False
        if not needs_np and not needs_pd:
            continue

        # insertion point: after docstring (if any), else before first stmt
        first = func.body[0]
        if is_docstring(first):
            insert_idx = first.end_lineno  # 0-indexed line after docstring
        else:
            insert_idx = first.lineno - 1  # 0-indexed: before first stmt

        # indentation from the first statement's line
        first_line = lines[first.lineno - 1]
        indent = " " * (len(first_line) - len(first_line.lstrip()))

        content: list[str] = []
        if needs_np:
            content.append(f"{indent}import numpy as np\n")
        if needs_pd:
            content.append(f"{indent}import pandas as pd\n")

        prepend_insert(insert_idx, content)

    # ── assemble result ───────────────────────────────────────────────────
    result: list[str] = []
    for i, line in enumerate(lines):
        if i in inserts_before:
            result.extend(inserts_before[i])
        if i not in lines_to_delete:
            result.append(line)
    if len(lines) in inserts_before:  # insertions after last line
        result.extend(inserts_before[len(lines)])

    new_source = "".join(result)
    if new_source == source:
        return None

    path.write_text(new_source, encoding="utf-8")
    return "CHANGED"


def main() -> None:
    src = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("src/seismicrna")

    changed = 0
    skipped: list[str] = []

    for path in sorted(src.rglob("*.py")):
        name = path.name
        if (
            name.endswith("_test.py")
            or name.endswith("_jit.py")
            or "__pycache__" in str(path)
        ):
            continue

        result = transform(path)
        rel = path.relative_to(src.parent.parent)
        if result == "CHANGED":
            print(f"  changed  {rel}")
            changed += 1
        elif result is not None:
            msg = result.replace("SKIP", "").strip()
            print(f"  skipped  {rel}  ({msg})")
            skipped.append(str(rel))

    print(f"\n{changed} files changed, {len(skipped)} skipped")
    if skipped:
        print("Files needing manual attention:")
        for p in skipped:
            print(f"  {p}")


if __name__ == "__main__":
    main()
