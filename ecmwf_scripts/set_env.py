import os
import subprocess

def set_environment():
    """
    This function sets the environment variables needed to execute the software.
    It assumes that you have the software installed somewhere in the $HOME directory.
    If not in $HOME, the program will ask you to provide a path where to look for it.
    Be sure to have permission to access this new path
    """

    path = [line for line in subprocess.check_output("find $HOME -iname 'molecularprofiles'",
                                                                          shell=True).splitlines()]
    if not path:
        new_path = raw_input('path to the molecularprofiles folder:')
        command = "find " + new_path +" -iname 'molecularprofiles'"
        path = [line for line in subprocess.check_output(command, shell=True).splitlines()]

    path_to_molecularprofiles = str(path[0])[2:-1]
    os.environ['MOLECULARPROFILES_DIR'] = path_to_molecularprofiles
    os.environ['ECMWF_DIR'] = os.environ['MOLECULARPROFILES_DIR'] + '/ecmwf_scripts/'
    os.environ['GDAS_DIR'] = os.environ['MOLECULARPROFILES_DIR'] + '/gdas_scripts/'
    os.environ['UTILS_DIR'] = os.environ['MOLECULARPROFILES_DIR'] + '/utils/'

    data_path = os.getcwd()
    os.environ['MOLECULAR_DATA_DIR'] = data_path

    print('Setting environment variables.')
    print('MOLECULARPROFILES_DIR: %s' % os.environ['MOLECULARPROFILES_DIR'])
    print('MOLECULAR_DATA_DIR: %s' % os.environ['MOLECULAR_DATA_DIR'])
    print('ECMWF_DIR: %s' % os.environ['ECMWF_DIR'])
    print('GDAS_DIR: %s' % os.environ['GDAS_DIR'])
    print('UTILS_DIR: %s' % os.environ['UTILS_DIR'])

if __name__ == '__main__':
    set_environment()