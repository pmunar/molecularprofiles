# molecularprofiles-init.sh sets some environment variables for molecularprofiles.
# 
# One way is to put this e.g. in ~/.bashrc or ~/.profile:
#
# export MOLEC_DIR=< location of your molecularprofiles software checkout >
# alias init_molec="source $MOLEC_DIR/molecularprofiles-init.sh"
#
# This way your python, PATH, PYTHONPATH, ... is not set to the
# molecularprofiles software when you log in to your system, yet you can
# conveniently set up your shell for molecularprofiles by executing the aliases:
# $ init_molec

echo "Adding molecularprofiles to PATH and PYTHONPATH"
export PATH=$PATH:$MOLEC_DIR/ecmwf_scripts:$MOLEC_DIR/gdas_scripts:$MOLEC_DIR/utils
export PYTHONPATH=$MOLEC_DIR/utils/:$MOLEC_DIR/ecmwf_scripts/:$MOLEC_DIR/gdas_scripts/:$PYTHONPATH

