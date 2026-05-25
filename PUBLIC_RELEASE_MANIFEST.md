# Public Release Manifest

Use `public_github_mumoas/` as the upload root for the journal-review GitHub package.

The public package should contain:

- `.gitignore`
- `README.md`
- `CODE_AVAILABILITY.md`
- `PUBLIC_RELEASE_MANIFEST.md`
- `LICENSE`
- `CITATION.cff`
- `pyproject.toml`
- `requirements.txt`
- `configs/`
- `examples/synthetic_data/.gitkeep`
- `mumoas/`
- `scripts/`
- `analysis_r/`
- `tests/`

The public package should not contain manuscript files, cover letters, disclosure forms, submitted figures, generated outputs, raw imaging data, masks, private clinical tables, omics matrices, single-cell matrices, Seurat objects, inferCNV objects, CellChat objects, model weights, checkpoints, private configuration files, or private source manifests.
