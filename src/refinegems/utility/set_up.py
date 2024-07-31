"""Collection of functions for setting up files and work environments.
"""

__author__ = 'Carolin Brune'

################################################################################
# requirements
################################################################################

import requests

from importlib.resources import files
from pathlib import Path
from tqdm import tqdm
from typing import Literal

################################################################################
# variables
################################################################################

PATH_MEDIA_CONFIG = files('refinegems.data.config').joinpath('media_config.yaml') #: :meta hide-value: 
# @TODO
# PATH_REFINEGEMS_CONFIG = files('data.config').joinpath('config.yaml')

################################################################################
# functions
################################################################################

# --------------------------
# download databases / files
# --------------------------
# @TEST
# @TODO : add an entry point?
def download_url(dowload_type:Literal['SwissProt gapfill'],
                 directory:str=None,k:int=10):
    
    # match URLS to type of database, that the user wants to download
    match dowload_type:
        case 'SwissProt gapfill':
            swissprot_api = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz'
            swissprot_mapping_api = 'https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cxref_brenda%2Cec%2Csequence&format=tsv&query=%28*%29+AND+%28reviewed%3Atrue%29'
            urls = {'SwissProt.fasta':swissprot_api, 'SwissProt_mapping.tsv':swissprot_mapping_api}
        case _:
            mes = f'Unknown database or file: {name}'
            raise ValueError(mes)

    # download each file
    for name,url in urls:
        r = requests.get(url, stream=True) # open download stream
        filename = Path(directory,name) if directory else Path(name)
        with open(filename, 'wb') as f:
            pbar = tqdm(desc=f'Downloading {name}', 
                        unit="B", unit_scale=True, unit_divisor=1024, 
                        total=int( r.headers['Content-Length'] )) # make the progress bar
            pbar.clear()  #  clear 0% info
            for chunk in r.iter_content(chunk_size=k*1024): 
                if chunk: # filter out keep-alive new chunks
                    pbar.update(len(chunk)) # update progress bar
                    f.write(chunk)
            pbar.close()

# ---------------------
# handling config files
# ---------------------

# @TODO : sth for gapfilling?
def download_config(filename:str='./my_config.yaml', type=Literal['media','refinegems']):
    """Load a configuration file from the package and save a copy of it for the user to edit.

    Args:
        - filename (str, optional): 
            Filename to write the config to/save it under as. 
            Defaults to './my_config.yaml'.
        - type (Literal['media','refinegems'], optional): 
            Type of configuration file to load.
            Can be 'media' for the media config file or 
            'refinegems' for the refinegmes pipeline configuration file. 
            Defaults to Literal['media','refinegems'].
    """

    def copy_config_yaml(infile:str, outfile:str):
        """Helper function for :py:func:`download_config`. 
        Performs the actual download of a yaml file into a copy for the user to edit.

        Args:
            - infile (str): 
                Path to the file to copy.
            - outfile (str): 
                Path to write the copy to.
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
        case 'gapfill':
            
            # copy_config_yaml(..., filename)
            pass

        # type not found
        case _:
            raise ValueError(F'Unknown type of config file detected: {type}')