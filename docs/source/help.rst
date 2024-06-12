Help and FAQ
============

This page provides help if you run into problems with ``refineGEMs``.

Help
----

Known Issues and Bugs
---------------------

.. hint:: 
    For bug reports please write issues on the `GitHub page <https://github.com/draeger-lab/refinegems/issues>`__ 
    or open a discussion `here <https://github.com/draeger-lab/refinegems/discussions>`__.

Pydantic
^^^^^^^^

Pydantic warning `underscore_attrs_are_private has been removed` has not - yet - caused any issues.
However, the core of the problem (= what causes the warning) has yet to be identified. 

FAQ
---

This site addresses frequently asked questions not addressed in other pages of the ``refineGEMs`` documentation.

How do I install ``refineGEMs``?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please see :ref:`Installation`.

How do I install ``refineGEMs`` as developer?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please see :ref:`Installation for developers`.

How do I cite ``refineGEMs``?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please see :ref:`How to cite`.

After restarting my device (MacOS) my conda set-up has the wrong python version. What can I do?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1. Deactivate all environments including the base environment
2. Activate base
3. Activate the ``refineGEMs`` environment

I cannot access ``python`` from within my conda environment in my VSCode terminal. What can I do?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Deactivate base and reactivate again:

.. code:: console

   conda deactivate
   conda deactivate
   conda activate base
   conda activate <your conda env>

.. note::
    TODO: Merge both FAQs for conda environment in VSCode and wrong Python version with MacOS?

My ``pipenv`` is not locking after f.ex. moving the repository. What can I do?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Try uninstalling ``pipenv`` and reinstalling it via pip. Then  run ``pipenv install`` and it should work again.

How to solve errors caused by ``pandoc`` (Development)?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
| `Here <https://stackoverflow.com/a/71585691>`__ is an answer that might help.
| If you install pandoc with conda use:

.. code:: console
    :class: copyable
    
    conda install -c conda-forge pandoc

How to solve errors caused by ``jinja2`` (Development)?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Just switch to version 3.0.3:

.. code:: console
    :class: copyable
    
    pip install jinja2==3.0.3

When do I use ``lab_strain=True``?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the parameter ``lab_strain`` is set to ``True`` it is assumed that the input model was created for a strain of an 
organism for which so far no data is available in any database. To get sufficient information for the knowledge-base 
stored within the model of the strain the ``PGAB`` pipeline described in ``SPECIMEN`` should be used before model 
creation. The resulting FASTA from the ``PGAB`` pipeline can then be used als input for CarveMe and is then also 
required as additional input for :py:mod:`~refinegems.curation.polish`. In this case the parameter ``lab_strain`` should 
be set to ``True``.
