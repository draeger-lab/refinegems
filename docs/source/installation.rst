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

``refineGEMs`` depends on the tools `MCC <https://github.com/Biomathsys/MassChargeCuration>`__ and 
`BOFdat <https://github.com/draeger-lab/BOFdat>`__ which cannot directly be installed via 
`PyPI <https://pypi.org/project/refineGEMs/>`__. Please install both tools before using ``refineGEMs``:

.. hint:: 

   BOFdat should be installed as stated below. The fork contains certain bug fixes such that the program is able 
   to run. If you have a functional version of BOFdat installed you can try to use it. In case you encounter problems 
   with your own version, please, install BOFdat as stated below.

.. code:: console
   :class: copyable

   pip install "masschargecuration@git+https://github.com/Biomathsys/MassChargeCuration@installation-fix"

.. code:: console
   :class: copyable

   pip install "bofdat@git+https://github.com/draeger-lab/BOFdat"


.. hint:: 

   For the installation for developers, refer to :ref:`installation for developers`