Installation
============

Installation via pip
--------------------
To install refinegmes as python package, simply install it via ``pip``:

.. code:: bash

   pip install refinegems

Installation via github
-----------------------

``refinegems`` is distributed via github repository, all
dependencies are denoted in a
`pipenv <https://pipenv.pypa.io/en/latest/>`__. You will need to install
``pipenv`` first. To install ``refinegems`` locally complete the
following steps:

.. code:: bash

   # install pipenv using pip
   pip install pipenv

   # clone or pull the latest source code
   git clone https://github.com/draeger-lab/gem_curation_template.git
   cd gem_curation_template

   # install all dependencies
   pipenv install

   # initiate a session in the virtual environment
   pipenv shell

The ``pipenv`` package can also be installed via Anaconda (recommended
if you are a Windows user).

Further Dependencies
--------------------

If you want to use the SBO terms annotation part you will need to
install PostgreSQL to your machine.

Make sure that you will be able to use psql from the command line.

The database containing the BiGG ID and EC number mappings

::

   sbo/create_dbs.sql

must be imported to a local PostgreSQL database to a selected user.

You can use the following command, run it from this directory:

.. code:: bash

   psql -U {your postgres username} -h localhost -d {your database name} < sbo/create_dbs.sql 

If you are a Windows user you will need to use a different command:
Enter into the psql shell by typing ``psql``, then create the database
with

.. code:: sql

   CREATE DATABASE sbo_ann;

Afterwards load the database with

.. code:: bash

   psql.exe -U postgres -d sbo_ann -f sbo\create_dbs.sql


Troubleshooting
---------------

-  If you get ``ImportError: DLL load failed while importing _sqlite3``
   when running main.py. Locate the ``sqlite3.dll`` file on you machine
   and add it to PATH.

-  If you use python 3.8 it everything should work, just edit the
   ``Pipfile`` entry to ``python_version = "3.8"`` before running
   ``pipenv install``.

-  If you can’t use ``psql``\ from the command line, a common issue is
   that its not added to PATH:

.. code:: bash

   locate psql | grep /bin
   export PATH={Output from the line above with /bin as line end}:$PATH

-  If you are a Windows user you will want to locate the installation
   manually. The path should look like something like this
   ``C:\Program Files\PostgreSQL\13\lib``.

-  If you run into
   ``psycopg2.OperationalError: fe_sendauth: no password supplied``:
   Change ``scram-sha256`` to ``trust`` in your file ``pg_hba.conf``
   (located probably in ``C:\Program Files\PostgreSQL\13\data``)

- If you run into a problem with ``pipenv`` not locking after f.ex. moving the repository try uninstalling ``pipenv`` and reinstalling it via pip. Then  run ``pipenv install`` and it should work again.


