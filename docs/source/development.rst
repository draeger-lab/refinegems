Development
===========

Development installation
------------------------

.. attention::
    The following packages need to be installed to be able to add content to the refineGEMs documentation.
    
    * `sphinx`
    * `nbsphinx`
    * `sphinx_rtd_theme`
    * `pandoc`
    * `ipython`

You can install the packages via pip to your local environment:

.. code:: bash

    pip install sphinx nbsphinx sphinx_rtd_theme pandoc ipython

If you run into an error with jinja2, just switch to version 3.0.3:

.. code:: bash
    
    pip install jinja2==3.0.3

Development notes
-----------------

You can enable debug logging by replacing ``level=logging.INFO``  with ``level=logging.DEBUG``.

If you want your print message to show in the log file, replace the ```print()`` statement by ``logging.info()``.