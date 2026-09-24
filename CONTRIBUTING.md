# Contributing

Open an issue or a pull request on this repository.

Run the tests before you send a change:

```bash
pip install -e .
PYTHONPATH=src python -m unittest tests.test_kernel
```

The license is Apache-2.0. Do not add a number the README does not name, and do not commit a private key or an environment file.

The pre-commit hook is `.githooks/pre-commit`. It refuses a private key, a token, or an env file, then runs the tests. Enable it once in a clone:

```bash
git config core.hooksPath .githooks
```
