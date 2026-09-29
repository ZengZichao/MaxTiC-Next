#!/bin/bash
# Bioconda build script for MaxTiC-Next (noarch Python package).
# This script is invoked by conda-build inside the build environment.
set -ex

# Install the package using pip (noarch Python, no compilation needed).
$PYTHON -m pip install . -vv --no-deps --no-build-isolation
