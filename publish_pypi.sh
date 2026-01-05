#!/bin/bash
set -e

# Cleanup previous builds
echo "Cleaning up previous builds..."
rm -rf dist/ build/ *.egg-info python/*.egg-info

# Install build tools if missing (assuming user has pip)
echo "Installing build tools..."
pip install --upgrade build twine

# Build
echo "Building package..."
python3 -m build

# Upload (using twine)
echo "Uploading to PyPI..."
# Note: User needs to have ~/.pypirc configured or provide credentials interactively
if [ -z "$(ls -A dist)" ]; then
   echo "Error: dist directory is empty!"
   exit 1
fi

twine upload dist/*
