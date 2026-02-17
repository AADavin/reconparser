# PyPI Release Guide

This document describes how to build and upload reconparser to PyPI.

## Prerequisites

1. Install build tools:
```bash
pip install build twine
```

2. Create accounts:
   - PyPI: https://pypi.org/account/register/
   - TestPyPI (for testing): https://test.pypi.org/account/register/

3. Configure API tokens:
   - Go to Account Settings → API tokens
   - Create a token for uploads
   - Store in `~/.pypirc`:

```ini
[pypi]
  username = __token__
  password = pypi-AgEIcHl...  # Your token here

[testpypi]
  username = __token__
  password = pypi-AgENdGVz...  # Your token here
```

## Pre-Release Checklist

1. **Update version** in `pyproject.toml`
2. **Update CHANGELOG** (if you create one)
3. **Run tests**:
   ```bash
   pytest tests/
   ```
4. **Check code formatting**:
   ```bash
   black src/
   isort src/
   flake8 src/
   ```
5. **Build documentation** (if applicable)
6. **Review README.md** for accuracy

## Building the Package

Clean previous builds:
```bash
rm -rf dist/ build/ *.egg-info
```

Build the package:
```bash
python -m build
```

This creates:
- `dist/reconparser-X.Y.Z.tar.gz` (source distribution)
- `dist/reconparser-X.Y.Z-py3-none-any.whl` (wheel)

## Testing on TestPyPI

Upload to TestPyPI first:
```bash
python -m twine upload --repository testpypi dist/*
```

Test installation:
```bash
pip install --index-url https://test.pypi.org/simple/ --no-deps reconparser
```

Test the installed package:
```python
from reconparser import ALEParser
parser = ALEParser("test.ale")
```

## Uploading to PyPI

Once testing is successful:

```bash
python -m twine upload dist/*
```

Verify on PyPI: https://pypi.org/project/reconparser/

Test installation:
```bash
pip install reconparser
```

## Version Numbering

Follow semantic versioning (MAJOR.MINOR.PATCH):

- **MAJOR**: Incompatible API changes
- **MINOR**: Add functionality (backwards-compatible)
- **PATCH**: Bug fixes (backwards-compatible)

Examples:
- `0.1.0` - Initial release
- `0.1.1` - Bug fix
- `0.2.0` - New features
- `1.0.0` - First stable release

## Post-Release

1. Create a git tag:
   ```bash
   git tag -a v0.1.0 -m "Release version 0.1.0"
   git push origin v0.1.0
   ```

2. Create a GitHub release with release notes

3. Update documentation if hosted separately

4. Announce on relevant platforms:
   - Twitter/X
   - Mastodon
   - Relevant mailing lists
   - Lab website/blog

## Troubleshooting

### Build fails
- Check `pyproject.toml` syntax
- Ensure all required files exist
- Verify package structure

### Upload fails
- Check API token
- Verify version not already uploaded
- Check internet connection

### Import fails after install
- Verify package structure
- Check dependencies installed
- Test in clean virtual environment

## Useful Commands

Check package metadata:
```bash
python -m twine check dist/*
```

List package contents:
```bash
tar tzf dist/reconparser-*.tar.gz
```

Show installed package info:
```bash
pip show reconparser
```
