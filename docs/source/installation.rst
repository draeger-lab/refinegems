Installation
============

The package ``refineGEMs`` and all its dependencies can be installed either via pip into virtual environments like 
``conda`` or ``pipenv`` or via Docker.
The latest version is tested with Python 3.10, 3.11 and 3.12.

.. hint::

   * For help and information about known bugs, refer to :ref:`Help and FAQ`.
   * For the installation for developers, refer to :ref:`installation for developers`

:icon:`fa-brands fa-python` Via pip
-----------------------------------

.. warning::

   | With the release of version 2.0.0, ``refineGEMs`` will only work with Python 3.10+.
   | For older versions Python 3.9 can still be used.

To install ``refineGEMs`` as Python package from `PyPI <https://pypi.org/project/refineGEMs/>`__, simply install it via ``pip``:

.. code:: console
   :class: copyable

   pip install refinegems

Some ``refineGEMs`` features require optional dependencies that are not needed
for the base installation:

.. code:: console
   :class: copyable

   pip install "refinegems[chebi]"

.. code:: console
   :class: copyable

   pip install "refinegems[ols]"

.. code:: console
   :class: copyable

   pip install "refinegems[sbo]"

To install all optional dependencies at once, run:

.. code:: console
   :class: copyable

   pip install "refinegems[optional]"

The connected tools `ModelPolisher <https://github.com/draeger-lab/ModelPolisher>`__\ :footcite:p:`Roemer2016_ModelPolisher,King2016_ModelPolisher`,
`MCC <https://github.com/draeger-lab/MassChargeCuration>`__\ :footcite:p:`Mostolizadeh2026_mcc` and
`BOFdat <https://github.com/draeger-lab/BOFdat>`__\ :footcite:p:`Lachance2019_bofdat` are also optional workflow
steps and currently need to be installed directly from GitHub before using the
corresponding functionality. If they are missing, ``refineGEMs`` reports the
missing dependency and skips the affected optional step where possible.

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


:icon:`fa-brands fa-docker` Via Docker
--------------------------------------

``refineGEMs`` can also be used via Docker. 
You can pull the latest image from (a) Docker Hub or (b) build it locally.

(a) Image from Docker Hub
^^^^^^^^^^^^^^^^^^^^^^^^^
To pull the image from Docker Hub, simply use:

.. code:: console
   :class: copyable

   docker pull biodatalab/refinegems:<tag>

(b) Local build
^^^^^^^^^^^^^^^
To build the Docker image locally, firstly clone the repository:

.. code:: console
   :class: copyable

   git clone "https://github.com/draeger-lab/refinegems.git"

Then change into the directory and build the image:

.. code:: console
   :class: copyable

   cd refinegems
   docker build -t refinegems .

The default image installs the runtime optional dependency group from
``pyproject.toml``, but excludes the documentation dependencies. Optional
connected tools that are currently installed directly from GitHub are included
by default and can be disabled for a smaller image:

.. code:: console
   :class: copyable

   # build without the optional connected GitHub tools
   docker build \
      --build-arg INSTALL_EXTERNAL_TOOLS=false \
      -t refinegems:runtime .

The full default can also be made explicit:

.. code:: console
   :class: copyable

   docker build \
      --build-arg INSTALL_EXTERNAL_TOOLS=true \
      -t refinegems:full .

How to use
^^^^^^^^^^

.. note::

   To provide the input files and retrieve the output files mount one folder as workspace folder to the Docker image 
   with ``-v``.

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

.. footbibliography::
