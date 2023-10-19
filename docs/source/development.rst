Development
===========

Development installation
------------------------

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
    

You can install the packages via pip to your local environment:

.. code:: bash

    pip install accessible-pygments sphinx nbsphinx sphinx_rtd_theme pandoc ipython sphinxcontrib-bibtex sphinx_copybutton

If you run into an error with jinja2, just switch to version 3.0.3:

.. code:: bash
    
    pip install jinja2==3.0.3

Development notes
-----------------

You can enable debug logging by replacing ``level=logging.INFO``  with ``level=logging.DEBUG``.

If you want your print message to show in the log file, replace the ```print()`` statement by ``logging.info()``.

Documentation notes
-------------------

We use the autoDocstring extension (njpwerner.autodocstring) for VSCode with the google format to generate function docstrings. To ensure a nice looking sphinx documentation, we add ``-`` to all variables that are passed as Args. And tuple returns are written as follows:

.. code:: python
    :linenos:

    """Description of the function...

    Args:
        - input1 (type): this is what input1 does

    Returns:
        tuple: Two tables (1) & (2)
            (1) pd.DataFrame: Table with charge mismatches
            (2) pd.DataFrame: Table with formula mismatches
    """

We are also trying to make input and return types explicit by declaring those in the function header:

.. code:: python
    :linenos:

    def my_func(input1: int, input2: str, input3: Model) -> tuple[str, int]:

More details for certain specifics can also be found `here <https://github.com/draeger-lab/refinegems/issues/74>`__.