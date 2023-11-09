Installation
============

Installation via pip
--------------------
To install refineGEMs as Python package from `PyPI <https://pypi.org/project/refineGEMs/>`__, simply install it via ``pip``:

.. code:: console
   :class: copyable

   pip install refineGEMs

The corresponding project site can be found `here <https://pypi.org/project/refineGEMs/>`__.

Installation via github
-----------------------

``refineGEMs`` is distributed via the GitHub repository.

The package and all its dependencies can be installed within virtual environments.

**Condaenv**

.. hint::
   | Installing refineGEMs in a conda environment on a  MacOS can cause the following issue:
   | -> After a restart of the device the whole conda set-up may have the wrong python version. 
   | If that is the case: 

   1. Deactivate all environments including the base environment
   2. Activate base
   3. Activate the refineGEMs environment

Setup a conda virtual environment and use its pip to install refineGEMs into that environment:

.. code:: bash

   # clone or pull the latest source code
   git clone https://github.com/draeger-lab/refinegems.git
   cd refinegems

   conda create -n <EnvName> python=3.9

   conda activate <EnvName>

   # check that pip comes from <EnvName>
   which pip

   pip install .

This will install all packags denoted in `setup.py`. 

If `which pip` does not show pip in the conda environment you can also create a local environment for which you can control the path and use its pip:

.. code:: bash

   conda create --prefix ./<EnvName>

   conda activate <path to EnvName>

   <EnvName>/bin/pip install .

**Pipenv**

You can use
`pipenv <https://pipenv.pypa.io/en/latest/>`__ to keep all dependencies together. You will need to install
``pipenv`` first. To install ``refineGEMs`` locally complete the
following steps:

.. code:: bash

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


Troubleshooting
---------------

-  If you get ``ImportError: DLL load failed while importing _sqlite3``
   when running main.py. Locate the ``sqlite3.dll`` file on you machine
   and add it to PATH.
- If you run into a problem with ``pipenv`` not locking after f.ex. moving the repository try uninstalling ``pipenv`` and reinstalling it via pip. Then  run ``pipenv install`` and it should work again.
- If you use VSCode terminals and have trouble accessing the python from within your conda environment, deactivate base and reactivate again:

.. code:: bash

   conda deactivate
   conda deactivate
   conda activate base
   conda activate <your conda env>
