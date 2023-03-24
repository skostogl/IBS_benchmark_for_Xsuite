#!/bin/bash

# Define a function to setup the package
function setup {
  if [ "$1" = "forMac" ]; then
    # Download and install miniconda
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-MacOSX-x86_64.sh -b -p ./miniconda -f
  elif [ "$1" = "forLinux" ]; then
    # Download and install miniconda
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    bash Miniconda3-latest-Linux-x86_64.sh -b -p ./miniconda -f
  fi

  if [ "$1" ]; then
    # Activate miniconda environment and create virtual environment
    source miniconda/bin/activate
  fi
  which python
  python -m pip install jupyterlab matplotlib pandas scipy numpy
  # Clone repository and install required packages
  git clone https://github.com/xsuite/tree_maker.git
  cd tree_maker
  git checkout release/v0.1.0
  python -m pip install .

}

# Call the function without arguments to setup the package
setup
# Call the function with an argument to download Python for a specific platform
#setup forMac
#setup forLinux
