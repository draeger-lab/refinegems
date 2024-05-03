Development
===========

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

.. code:: bash
    :class: copyable

    pip install accessible-pygments sphinx nbsphinx sphinx_rtd_theme pandoc ipython sphinxcontrib-bibtex sphinx_copybutton

.. code:: bash
    :class: copyable

    python -m pip install pip-tools

Troubleshooting of installation issues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| If you run into an error with pandoc, `here <https://stackoverflow.com/a/71585691>`__ is an answer that might help.
| If you install pandoc with conda use:

.. code:: bash
    :class: copyable
    
    conda install -c conda-forge pandoc

If you run into an error with jinja2, just switch to version 3.0.3:

.. code:: bash
    :class: copyable
    
    pip install jinja2==3.0.3

Updating the `requirements.txt`
-------------------------------
| To create the `requirements.txt` adjust the `requirements.in` file as needed in the folder docs.
| Then navigate to the folder docs in the command line:

.. code-block:: bash
    :class: copyable

    cd docs

and use the following command to automatically generate the new `requirements.txt`:

.. code-block:: bash
    :class: copyable
    
    python3.9 -m piptools compile --output-file=requirements.txt requirements.in

Debugging switches
------------------

- You can enable debug logging by replacing ``level=logging.INFO``  with ``level=logging.DEBUG``.
- If you want your print message to show in the log file, replace the ```print()`` statement by ``logging.info()``.
- For debugging of pandas warnings or issues ``pd.options.mode.chained_assignment = None`` needs to be commented out.

Guidelines for code documentation
---------------------------------

We use the autoDocstring extension (njpwerner.autodocstring) for VSCode with the google format to generate function docstrings. 
To ensure a nice looking sphinx documentation, we add ``-`` to all variables that are passed as Args. And tuple returns are written as follows:

If you use VSCode, a mustache file for the documentation style that can be integrated into VSCode can be found in the ``dev``` directory.

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