Installation
============

The package ``refineGEMs`` and all its dependencies can be installed within virtual environments like ``conda`` or ``pipenv``.

.. warning::

   | With the release of version 2.0.0, refineGEMs will only work with Python 3.10+.
   | For older versions Python 3.9 can still be used.

.. hint::

   For help and information about known bugs, refer to :ref:`Help and FAQ`.

To install refineGEMs as Python package from `PyPI <https://pypi.org/project/refineGEMs/>`__, simply install it via ``pip``:

.. code:: console
   :class: copyable

   pip install refineGEMs

The corresponding project site can be found `here <https://pypi.org/project/refineGEMs/>`__.

``refineGEMs`` depends on the tools `ModelPolisher <https://github.com/draeger-lab/ModelPolisher>`__, `MCC <https://github.com/Biomathsys/MassChargeCuration>`__ and 
`BOFdat <https://github.com/draeger-lab/BOFdat>`__ which cannot directly be installed via 
`PyPI <https://pypi.org/project/refineGEMs/>`__. Please install all three tools before using ``refineGEMs``:

.. hint:: 

   BOFdat should be installed as stated below. The fork contains certain bug fixes such that the program is able 
   to run. If you have a functional version of BOFdat installed you can try to use it. In case you encounter problems 
   with your own version, please, install BOFdat as stated below.

.. code:: console
   :class: copyable

   pip install "model-polisher@git+https://github.com/draeger-lab/MPClient"

.. code:: console
   :class: copyable

   pip install "masschargecuration@git+https://github.com/draeger-lab/MassChargeCuration"

.. code:: console
   :class: copyable

   pip install "bofdat@git+https://github.com/draeger-lab/BOFdat"

``refineGEMs`` can also be used via Docker. To build the docker image, firtsly clone the repository:

.. code:: console
   :class: copyable

   git clone "https://github.com/draeger-lab/refinegems.git"

Then change into the directory and build the image:

.. code:: console
   :class: copyable

   cd refinegems \
   docker build -t refinegems .

.. note:: 

   To provide the input files and retrieve the output files mount one folder as workspace folder to the Docker image 
   with `-v`.

The default command executed by the image is ``refinegems -h`` and provides the help information for the CLI of 
``refineGEMs``.

.. code:: console
   :class: copyable

   docker run refinegems -h

To use the image interactively and open a bash shell, run the following command:

.. code:: console
   :class: copyable

   docker run -it --entrypoint bash refinegems

To use the image for specific commands, you can simply use every of the CLI commands as entrypoint. 
For example, to curate a (draft) model, run:

.. code:: console
   :class: copyable

   docker run --name <container_name> -v <user_folder>:/rg_cont refinegems analyse stats ./path/to/model.xml


.. hint::

   For the installation for developers, refer to :ref:`installation for developers`