"""Collection of functions for setting up files and work environments."""

__author__ = "Carolin Brune"

################################################################################
# requirements
################################################################################

import requests
import subprocess
import warnings

from importlib.resources import files
from pathlib import Path
from tqdm import tqdm
from typing import Literal

################################################################################
# variables
################################################################################

PATH_MEDIA_CONFIG = files("refinegems.data.config").joinpath(
    "media_config.yml"
)  #: :meta hide-value:

################################################################################
# functions
################################################################################

# --------------------------
# download databases / files
# --------------------------


def download_url(
    download_type: Literal["SwissProt gapfill"],
    directory: str = None,
    k: int = 10,
    t: int = 1,
):
    """Download files necessary for certain functionalities of the toolbox from
    the internet.

    Currently available:

    - 'SwissProt gapfill': download files needed for the :py:class:`~refinegems.classes.gapfill.GeneGapFiller`


    Args:
        - dowload_type (Literal['SwissProt gapfill']):
            Type of files to download.
        - directory (str, optional):
            Path to a directory to save the downloaded files to.
            Defaults to None (So the current working directory is used).
        - k (int, optional):
            Chunksize in kB.
            Defaults to 10.
        - t (int, optional):
            Number of threads to use for some additional setups, e.g.
            DIAMOND database creation.
            Defaults to 1.

    Raises:
        - ValueError: Unknown database or file
    """

    # Ensure that directory is not None
    if directory is None:
        directory = Path.cwd()

    # match URLS to type of database, that the user wants to download
    match download_type:
        case "SwissProt gapfill":
            swissprot_api = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
            swissprot_mapping_api = "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cxref_brenda%2Cec&format=tsv&query=%28*%29+AND+%28reviewed%3Atrue%29"
            urls = {
                "SwissProt.fasta": swissprot_api,
                "SwissProt_mapping.tsv": swissprot_mapping_api,
            }
            #   1.: TSV with UniprotID, BRENDA and EC -7.7MB (26.07.2024)
            #   2.: FASTA with sequences ~280MB (26.07.2024)
        case _:
            mes = f"Unknown database or file: {download_type}"
            raise ValueError(mes)

    # download each file
    for name, url in urls.items():
        r = requests.get(url, stream=True)
        filename = Path(directory, name)

        # Check if Content-Length is available
        total_length = r.headers.get("Content-Length")  # Make the progress bar
        if total_length is None:
            # Content-Length is missing, so we download without a progress bar
            with open(filename, "wb") as f:
                for chunk in r.iter_content(chunk_size=k * 1024):
                    if chunk:
                        f.write(chunk)
        else:
            total_length = int(total_length)
            with open(filename, "wb") as f:
                pbar = tqdm(
                    desc=f"Downloading {name}",
                    unit="B",
                    unit_scale=True,
                    unit_divisor=1024,
                    total=total_length,
                )
                pbar.clear()
                for chunk in r.iter_content(chunk_size=k * 1024):
                    if chunk:  # Filter out keep-alive new chunks
                        pbar.update(len(chunk))  # Update progress bar
                        f.write(chunk)
                pbar.close()

    # additional setups
    match download_type:
        # SwissProt gapfill
        case "SwissProt gapfill":
            # Create db folder for DIAMONd if non-existent
            try:
                Path(directory, "db").mkdir(parents=True, exist_ok=False)
                print("Creating new directory " + str(Path(directory, "db")))
            except FileExistsError:
                warnings.warn(
                    "Given directory already exists. High possibility of files being overwritten."
                )
            Path(directory, "db").mkdir(parents=True, exist_ok=True)
            # create DIAMOND database
            print(f"create DIAMOND database for {download_type} using:")
            print(
                f'diamond makedb --in {Path(directory, "SwissProt.fasta")} --db {str(Path(directory,"db","SwissProt.dmnd"))} --threads {int(t)}'
            )
            subprocess.run(
                [
                    "diamond",
                    "makedb",
                    "--in",
                    str(Path(directory, "SwissProt.fasta")),
                    "--db",
                    str(Path(directory, "db", "SwissProt.dmnd")),
                    "--threads",
                    str(t),
                ]
            )
        # Type for which no extra setup is needed
        case _:
            pass


# ---------------------
# handling config files
# ---------------------


def download_config(filename: str = "./my_config.yaml", type=Literal["media"]):
    """Load a configuration file from the package and save a copy of it for the user to edit.

    Args:
        - filename (str, optional):
            Filename to write the config to/save it under as.
            Defaults to './my_config.yaml'.
        - type (Literal['media'], optional):
            Type of configuration file to load.
            Can be 'media' for the media config file.
            Defaults to Literal['media'].
    """

    def copy_config_yaml(infile: str, outfile: str):
        """Helper function for :py:func:`download_config`.
        Performs the actual download of a yaml file into a copy for the user to edit.

        Args:
            - infile (str):
                Path to the file to copy.
            - outfile (str):
                Path to write the copy to.
        """
        with open(infile, "r") as cfg_file, open(outfile, "w") as cfg_out:
            for line in cfg_file:
                cfg_out.write(line)

    # copy an examplary version of the config file for the user to edit it
    match type:

        # copy media config
        case "media":
            copy_config_yaml(PATH_MEDIA_CONFIG, filename)

        # type not found
        case _:
            raise ValueError(f"Unknown type of config file detected: {type}")
