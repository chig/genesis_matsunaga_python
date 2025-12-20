#!/bin/bash

# Change to the script's directory (python_interface/)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Add parent directory to PYTHONPATH so python can find the python_interface module
export PYTHONPATH="${SCRIPT_DIR}/..:${PYTHONPATH}"

python -m python_interface.test_crd_convert "$@"
python -m python_interface.test_trj "$@"
python -m python_interface.test_wham "$@"
python -m python_interface.test_mbar_1d "$@"
python -m python_interface.test_mbar_block "$@"
python -m python_interface.test_avecrd "$@"
python -m python_interface.test_kmeans "$@"
python -m python_interface.test_hb_atom "$@"
python -m python_interface.test_hb_snap "$@"
python -m python_interface.test_rmsd "$@"
python -m python_interface.test_drms "$@"
python -m python_interface.test_rg "$@"
python -m python_interface.test_msd "$@"
python -m python_interface.test_diffusion "$@"
python -m python_interface.test_mdanalysis "$@"
python -m python_interface.test_mdtraj "$@"
