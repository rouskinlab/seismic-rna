# Instructions for AI coding agents

## Releasing a new version

Before doing ANY work related to releasing SEISMIC-RNA — on GitHub, PyPI, or
Bioconda; including drafting release notes, bumping the version, editing
`meta.yaml`, or diagnosing a failed release workflow — read
[`RELEASE_CHECKLIST.md`](RELEASE_CHECKLIST.md) in full first and follow it.
Releasing is a multi-platform process with steps that are easy to silently
skip (e.g. Bioconda's `meta.yaml` does not automatically stay in sync with
`pyproject.toml`), so do not rely on memory or partial context from the
current conversation — always re-read the checklist file itself, since it is
kept up to date independently of any one conversation.
