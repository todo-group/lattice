# Publishing to PyPI

The PyPI distribution name is `lattice-graph-core`; the Python import name is
`lattice`.

## Local release check

```sh
python3 -m venv .venv
.venv/bin/python -m pip install maturin numpy twine
.venv/bin/python -m maturin build --release --out dist
.venv/bin/python -m maturin sdist --out dist
.venv/bin/python -m twine check dist/*
.venv/bin/python -m pip install --force-reinstall dist/lattice_graph_core-*.tar.gz
```

## Upload

Preferred: create a GitHub release or manually run the `python-publish`
workflow. Configure a PyPI Trusted Publisher for:

* repository: `todo-group/lattice`
* workflow: `python-publish.yml`
* environment: `pypi`
* project: `lattice-graph-core`

Manual upload is also possible:

```sh
.venv/bin/python -m maturin upload dist/*
```
