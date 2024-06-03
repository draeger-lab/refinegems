Usage
======

This page show to to quickly access the ``refinegems`` toolbox. More more detailed examples see :ref:`About different applications`.

Access points
-------------

via command line/terminal
^^^^^^^^^^^^^^^^^^^^^^^^^

The main functionalities of refineGEMs can be accessed directly from the command line after a successfull installation.

For more information, call the following in your command line or terminal

.. code-block:: bash

  refinegems --help

or refer to the following section: :ref:`Command line access points`.

access from inside Python
^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, refineGEMs can be directly imported into Python or a python script for full access of all modules, classes and more.

.. code-block:: python 

  import refinegems as rg

Or, e.g.:

.. code-block:: python

  from refinegems.analysis.growth import growth_analysis
  

The media database  
------------------

The media database can be accessed via a mutitude of different routes.

.. hint::

  For more information about the database and its contents, refer to :ref:`Media & Subsets`.

Direct access via SQL
^^^^^^^^^^^^^^^^^^^^^

One way to inspect the database is to download it from the `GitHub page <https://github.com/draeger-lab/refinegems>`__
and open it with an IDE with an SQL-extension (e.g. VSCode) or with a graphical SQL interface like SQLite.

Access via the command line
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``refineGEMs`` command line tool provides limited access to some of the basic functionalities concering
the media database, including downloading a config file and retrieving a list of all available media:

.. code-block:: bash

  refinegems media --help

For more information, see section :ref:`refinegems media`.

Access via Python
^^^^^^^^^^^^^^^^^

To access the media from inside a Python script, the :py:mod:`~refinegems.classes.medium`
provides classes and functions.

For example, to extract the medium `LB` run:

.. code-block:: python

  from refinegems.classes import Medium

  lb = medium.load_medium_from_db('LB')


Access the database via a config.yaml
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If the idea is to load a list of media - from the database or externel, with or without modifications -
it might be useful to use the media configuration file. There, different media names with options for 
modifications can be listed and subsequently loaded using :py:class:`~refinegems.classes.medium` or used as input for other 
entry points that requires a list of media or the config itself. 

.. toctree::

  The configuration file <config-desc/media-yaml>
  How to adjust the config <config-desc/media-yaml-howto>



