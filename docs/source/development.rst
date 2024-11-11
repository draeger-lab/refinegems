Development
===========

To maintain or extend the toolbox ``refineGEMs`` please install the package via GitHub.

.. warning::
   refineGEMs requires at least Python 3.10 since version 2.0.0.

.. hint::

   For help and information about known bugs, refer to :ref:`Help and FAQ`.

Installation for developers
---------------------------

``refineGEMs`` depends on the tools `MCC <https://github.com/Biomathsys/MassChargeCuration>`__ and 
`BOFdat <https://github.com/draeger-lab/BOFdat>`__ which cannot directly be installed via 
`PyPI <https://pypi.org/project/refineGEMs/>`__ or the `pyproject.toml`. Please install both tools before using 
``refineGEMs`` into the corresponding environment:

.. code:: console
    :class: copyable

    pip install "masschargecuration@git+https://github.com/Biomathsys/MassChargeCuration@installation-fix"

.. code:: console
    :class: copyable

    pip install "bofdat@git+https://github.com/draeger-lab/BOFdat"


Into a Conda environment
^^^^^^^^^^^^^^^^^^^^^^^^

Setup a conda virtual environment and use its pip to install ``refineGEMs`` into that environment.

.. code:: console

   # clone or pull the latest source code
   git clone https://github.com/draeger-lab/refinegems.git
   cd refinegems

   conda create -n <EnvName> python=<Specific Python version >= 3.10>

   conda activate <EnvName>

   # check that pip comes from <EnvName>
   which pip

   pip install .

This will install all packags denoted in `pyproject.toml`. 

If `which pip` does not show pip in the conda environment you can also create a local environment for which you can 
control the path and use its pip:

.. code:: console

   conda create --prefix ./<EnvName>

   conda activate <path to EnvName>

   <EnvName>/bin/pip install .

Into a Pipenv environment
^^^^^^^^^^^^^^^^^^^^^^^^^

You can use `pipenv <https://pipenv.pypa.io/en/latest/>`__ to keep all dependencies together. Therefore, you will need 
to install ``pipenv`` first. To install ``refineGEMs`` locally complete the following steps:

.. code:: console

   # install pipenv using pip
   pip install pipenv

   # clone or pull the latest source code
   git clone https://github.com/draeger-lab/refinegems.git
   cd refinegems

   # install all dependencies from Pipfile
   pipenv install .

   # initiate a session in the virtual environment
   pipenv shell

The ``pipenv`` package can also be installed via Anaconda (recommended
if you are a Windows user).

.. hint::

   If you want to be able to savely import the package from anywhere while also retaining the possibility to edit the 
   code, it is recommended to change the :code:`pip install` line from the code blocks to 
   :code:`pip install -e . --config-settings editable_mode=strict`.

Additional packages required for development
--------------------------------------------

.. attention::
    The following packages need to be installed to be able to add content to the refineGEMs documentation.
    
    * `accessible-pygments`
    * `ipython`
    * `nbsphinx`
    * `pandoc`
    * `sphinx`
    * `sphinx_copybutton`
    * `sphinx_rtd_theme`
    * `sphinxcontrib-bibtex`
    

    In addition, `pip-compile` should be installed to update the `requirements.txt` for the next release.

Installing the packages
^^^^^^^^^^^^^^^^^^^^^^^
You can install the packages via pip to your local environment:

.. code:: console
    :class: copyable

    pip install accessible-pygments sphinx nbsphinx sphinx_rtd_theme pandoc ipython sphinxcontrib-bibtex sphinx_copybutton

.. code:: console
    :class: copyable

    python -m pip install pip-tools

Alternatively, install the tool with the extra `docs`, e.g. 

.. code:: console
    :class: copyable

     pip install -e ".[docs]" --config-settings editable_mode=strict  

Updating the `requirements.txt`
-------------------------------
| To create the `requirements.txt` adjust the `requirements.in` file as needed in the folder docs.
| Then navigate to the folder docs in the command line:

.. code-block:: console
    :class: copyable

    cd docs

and use the following command to automatically generate the new `requirements.txt`:

.. code-block:: console
    :class: copyable
    
    python3 -m piptools compile --strip-extras --output-file=requirements.txt requirements.in

To bump to the newest versions possible, use the following command in the `docs` directory:

.. code-block:: console
    :class: copyable

    pip-compile --upgrade


Debugging switches
------------------

- You can enable debug logging by replacing ``level=logging.INFO``  with ``level=logging.DEBUG``.
- If you want your print message to show in the log file, replace the ```print()`` statement by ``logging.info()``.
- For debugging of pandas warnings or issues ``pd.options.mode.chained_assignment = None`` needs to be commented out.
- | Additionally, some modules contain comment blocks inf the format shown below.
  | By enabling the code lines between the dotted lines, a debugging-mode is run, which e.g. subsets the data to shorten the runtime to make debugging faster.

.. code-block:: python
    :class: copyable

    # @DEBUG ...............
    # some code
    # ......................

Guidelines for code documentation
---------------------------------

We use the autoDocstring extension (njpwerner.autodocstring) for VSCode with the google format to generate function docstrings. 
To ensure a nice looking sphinx documentation, we add ``-`` to all variables that are passed as Args. And tuple returns are written as follows:

If you use VSCode, a `mustache file for the documentation style <https://github.com/draeger-lab/refinegems/blob/dev-2/dev/docstring-format.mustache>`__ that can be integrated into VSCode (`dev <https://github.com/draeger-lab/refinegems/tree/dev-2/dev>`__ directory of refineGEMs).

.. code:: python
    :linenos:

    # Tuple output & Single input 
    """Description of the function...

    Args:
        - input1 (type):    
            this is what input1 does

    Returns:
        tuple: 
            Two tables (1) & (2)

            (1) pd.DataFrame: Table with charge mismatches
            (2) pd.DataFrame: Table with formula mismatches
    """

    # Single output with multiple possibilities & multiple inputs
    """Description of the function...

    Args:
        - input1 (type): 
            this is what input1 does
        - input2 (type):
            this is what input2 does
        - input3 (type): 
            this is what input3 does

    Returns:
        (1) Case: str

            Return value 1

        (2) Case: np.nan
        
            Return value 2
    """

We are also trying to make input and return types explicit by declaring those in the function header:

.. code:: python
    :linenos:

    def my_func(input1: int, input2: str, input3: Model) -> tuple[str, int]:

More details for certain specifics can also be found `here <https://github.com/draeger-lab/refinegems/issues/74>`__.


Information about working on the media database
-----------------------------------------------

Add of update information in the database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At the end of the :py:mod:`~refinegems.classes.medium` module are a set of funtions for automatic curation 
of the database. 

More information about how to run these can be found in the ``db_extension.ipynb`` notebook in the `dev <https://github.com/draeger-lab/refinegems/tree/dev-2/dev>`__ folder 
inside the `GitHub repository <https://github.com/draeger-lab/refinegems/>`_.

Create docs for media and subsets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After adding a new medium or subset to the database or updating existing information, 
the new documentation pages (.rst) can be generated automatically

For the media definition, use :py:func:`~refinegems.classes.medium.Medium.export_to_file`:

.. code:: python
    :linenos:

    new_medium = load_medium_from_db(<name>)
    new_medium.export_to_file(type='docs',dir=<path>)

To create the documentation page for a subset, use :py:func:`~refinegems.classes.medium.generate_docs_for_subset`:

.. code:: python 
    :linenos:

    generate_docs_for_subset(<name>,folder='<path>')