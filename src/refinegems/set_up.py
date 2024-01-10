"""Collection of functions for setting up files and work environments.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

from importlib.resources import files
# from tqdm import tqdm
from typing import Literal

################################################################################
# variables
################################################################################

PATH_MEDIA_CONFIG = files('data.config').joinpath('media_config.yaml')
# @TODO
# PATH_REFINEGEMS_CONFIG = files('data.config').joinpath('config.yaml')

################################################################################
# functions
################################################################################

# ---------------------
# handling config files
# ---------------------

def download_config(filename:str='./my_config.yaml', type=Literal['media','refinegems']):
    """Load a configuration file from the package and save a copy of it for the user to edit.

    Args:
        filename (str, optional): Filename to write the config to/save it under as. 
            Defaults to './my_config.yaml'.
        type (Literal['media','refinegems'], optional): Type of configuration file to load.
            Can be 'media' for the media config file or 
            'refinegems' for the refinegmes pipeline configuration file. 
            Defaults to Literal['media','refinegems'].
    """

    def copy_config_yaml(infile:str, outfile:str):
        """Helper function for :py:func:`download_config`. 
        Performa the actual download of a yaml file into a copy for the user to edit.

        Args:
            infile (str): Path to the file to copy.
            outfile (str): Path to write the copy to.
        """
        with open(infile, "r") as cfg_file, open(outfile, 'w') as cfg_out:
                for line in cfg_file:
                    cfg_out.write(line)

    # copy an examplary version of the config file for the user to edit it
    match type:

        # copy media config
        case 'media':           
            copy_config_yaml(PATH_MEDIA_CONFIG, filename)

        # @TODO
        # refinegems config
        case 'refinegems':
            
            # copy_config_yaml(PATH_REFINEGEMS_CONFIG, filename)
            pass

        # type not found
        case _:
            raise ValueError(F'Unknown type of config file detected: {type}')