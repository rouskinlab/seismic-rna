# Release Checklist

Follow this checklist, in order, every time SEISMIC-RNA is released. It is the
authoritative, up-to-date quick reference; `src/devdocs/tex/release.tex`
has the fuller narrative explanation of *why* each step exists, but this file
is what to actually follow — update this file (and release.tex, if the
process itself changes) whenever a step here turns out to be wrong or
outdated.

## 0. Before releasing anything

- [ ] Version bumped in `src/seismicrna/core/version.py` (`__version__`).
- [ ] Version bumped in `meson.build`. (Unfortunately, the version must also
      be defined this file and MUST match `version.py`.)
- [ ] `CHANGELOG.md` updated with the new version's changes.
- [ ] All GitHub Actions workflows on the latest commit to `main` are green
      (check the "Actions" tab — there is more than one workflow per commit,
      check ALL of them, not just the first).
- [ ] No uncommitted or unpushed changes you meant to include.

## 1. Release on GitHub

- [ ] Open the repo → "Releases" → "Create a new release".
- [ ] Tag: exactly `v` followed by the version number, e.g. `v0.27.0`. This
      exact format is required — Bioconda's autobump bot parses it, and a
      malformed tag silently breaks the Bioconda step below.
- [ ] Release title: short (≤10 words), summarizes the release.
- [ ] Release description: starts with "What's new in x.y.z", bullet points
      by category (features, fixes, breaking changes), then click "Generate
      release notes" at the bottom to append the auto-generated commit/PR
      summary.
- [ ] Publish. This triggers `.github/workflows/build-publish.yml`.

## 2. Release on PyPI

- [ ] Confirm `build-publish.yml` completes the sdist build and both wheel
      jobs (Linux, macOS) with green checkmarks. If a build fails, fix the
      root cause and re-run ("Re-run all jobs") rather than working around it.
- [ ] The final publish step requires manual approval (a repository owner —
      currently `matthewfallan` or `justinaruda`) in the `publish-pypi`
      GitHub Environment. Approve only after confirming tests passed and the
      code is functional.
- [ ] Once approved, PyPI publishing happens automatically via Trusted
      Publisher Management. Verify with `pip install seismic-rna==X.Y.Z` in a
      clean environment.
- [ ] If GitHub was already released but PyPI publishing needs to be
      re-triggered without creating a new GitHub release, re-run
      `build-publish.yml` manually via `workflow_dispatch` (Actions tab →
      "Build and Publish Project on PyPI" → "Run workflow"), selecting the
      release tag as the ref — NOT `main`, or you may publish the wrong code
      to a version number that's already taken.

## 3. Release on Bioconda

This step lags the other two and needs the most manual attention, because
there is currently **no automation that keeps `meta.yaml`'s dependency pins
in sync with `pyproject.toml`** (the old `make_conda_recipe.py` generator was
deleted from this repo; see `src/devdocs/tex/release.tex` for history). Do
not assume the recipe is correct just because a PR exists.

- [ ] **Before bumping any dependency floor in `pyproject.toml`, confirm it's
      actually required** (a bug fix, a new API you use), not just "the
      newest version available that day." A too-tight floor can pass PyPI
      (which always has the newest release) while silently breaking
      Bioconda, since conda-forge lags PyPI by anywhere from days to months.
      This has bitten two releases in a row: the `>=3.13.15` Python floor
      broke the Ubuntu wheel build, and `matplotlib`/`plotly`/
      `pyahocorasick`/`scipy` floors bumped in one "update everything to
      latest" commit broke this Bioconda build because conda-forge didn't
      yet have those exact versions (and `matplotlib>=3.11` additionally
      pulled in `libraqm`, which needs a newer glibc than bioconda's build
      image provides). When in doubt, relax the floor to whatever the code
      actually exercises rather than to the latest release.
- [ ] Within about a day of the GitHub tag going out, Bioconda's autobump bot
      should open a PR at
      https://github.com/bioconda/bioconda-recipes/pulls titled
      "Update seismic-rna to X.Y.Z" (branch `bump/seismic_rna`). If it hasn't
      appeared after a day or two, it can be created manually — see
      `src/devdocs/tex/release.tex` for the manual-fork procedure.
- [ ] **Do not assume the bot's PR is sufficient.** The bot only updates the
      `version` and the source tarball `sha256` in `recipes/seismic-rna/meta.yaml`.
      It does NOT update dependency versions, add new dependencies, remove
      dropped dependencies, or update the required Python version.
- [ ] Manually diff `recipes/seismic-rna/meta.yaml`'s `requirements` section
      against:
    - `pyproject.toml`'s `dependencies` list (Python packages) and
      `requires-python` (must match the `python` pin in `host`/`run`, and the
      `skip: True  # [py < NNN or py >= NNN]` build selector).
    - `environment.yml`'s non-Python system dependencies (bowtie2, fastp,
      rnastructure, samtools, seqkit, viennarna, etc.) and their version
      floors.
    - The actual `import`/subprocess dependencies in the code (e.g. grep for
      new external tools) if a new release added or removed one — `meta.yaml`
      dependencies are not purely a mirror of `pyproject.toml`, since it also
      lists non-Python runtime tools like `maven`, `jgo`, and `openjdk` that
      `pyproject.toml` doesn't cover.
    - Fix any mismatches by editing `meta.yaml` directly and pushing to the
      PR branch (if you have push access) or via a new PR against the bot's
      branch.
- [ ] Lint locally before pushing, if you have the `bioconda` Conda
      environment set up (`conda create -n bioconda -c conda-forge -c
      bioconda bioconda-utils`):
      `bioconda-utils lint --packages seismic-rna` and
      `bioconda-utils build --packages seismic-rna` should both report OK.
- [ ] Watch the PR's CI checks (Linux + macOS builds/tests) turn green. If a
      build fails, read the actual error in the job log — the summary badge
      alone is not enough; conda-build failures are often several hundred
      lines above the final traceback.
- [ ] Once CI is green, wait for a Bioconda maintainer to merge (can take
      several days). Check back periodically.
- [ ] After merge, verify with `conda install -c bioconda
      seismic-rna==X.Y.Z` (may take a few hours to propagate after merge).
- [ ] Delete the now-unneeded local/remote release branches.
