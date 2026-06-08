# Building the Documentation

This documentation is a [Jupyter Book](https://jupyterbook.org). You can build it
locally to preview changes before pushing. The only dependency is `jupyter-book`
itself (see `docs/requirements.txt`).

```{tip}
Copy each command on its own, without the surrounding prose. Some shells (e.g.
zsh on macOS) do **not** treat `#` as an inline comment and will pass it to the
command as an argument.
```

## One-time setup

In your local Python environment, install the docs dependency:

```bash
pip install -r docs/requirements.txt
```

(Equivalently: `pip install "jupyter-book>=0.15,<1.0"`.)

## Build

From the repository root:

```bash
jupyter-book build docs/
```

`jb build docs/` is a shorter alias for the same command.

## View

Open the generated HTML in a browser. On macOS:

```bash
open docs/_build/html/index.html
```

On Linux:

```bash
xdg-open docs/_build/html/index.html
```

Or paste the `file://` path to `docs/_build/html/index.html` into your browser bar.

## Tips

- **Clean rebuild.** The build caches; if something looks stale (especially after
  editing `_toc.yml`), force a full rebuild with `jupyter-book build docs/ --all`,
  or delete `docs/_build/` first.
- **`docs/_build/` is gitignored** — it's regenerated locally and never committed,
  so it won't appear when you pull.
- **No code runs** during the build (`execution_mode: off` in `_config.yml`), so it's
  fast and needs neither the model nor any data — it renders the Markdown and the
  notebook *source* only.
- **Hosted build.** The same book builds on ReadTheDocs via `.readthedocs.yml`; the
  local command above is all you need to preview your edits.

## Editing the book

- Pages are Markdown (`.md`) and notebooks (`.ipynb`) under `docs/`.
- The table of contents / page order is set in `docs/_toc.yml`.
- Book-wide settings (title, theme, etc.) are in `docs/_config.yml`.
- To add a page, create the file and add its path to `docs/_toc.yml`. Every page
  needs a top-level `#` title, or its TOC link won't be generated.
