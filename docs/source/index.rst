Welcome to refineGEMs!
======================================
``refineGEMs`` is a Python-based toolbox for the curation and analysis of genome-scale metabolic models (GEMS).

Overview
--------

This documentation provides information on how to install the toolbox, how to use it - either via 
command line or code - and provides different examples of use cases. Furthermore, it provides detailed 
information about the in-build media database and the code behind the modules.

Quick access to main sites
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: images/buttons/refineGEMs_button_install.png
  :height: 0px
  :width: 0px

.. raw:: html

  <a class="reference external image-reference" href="installation.html">
    <img src='_images/refineGEMs_button_install.png' alt='Installation' title='Installation' style="width: 20%;">
  </a>

.. image:: images/buttons/refineGEMs_button_usage.png
  :height: 0px
  :width: 0px

.. raw:: html

  <a class="reference external image-reference" href="main-modules.html">
    <img src='_images/refineGEMs_button_usage.png' alt='Applications' title='Applications' style="width: 20%;">
  </a>

.. image:: images/buttons/refineGEMs_button_media.png
  :height: 0px
  :width: 0px

.. raw:: html

  <a class="reference external image-reference" href="database/intro.html">
    <img src='_images/refineGEMs_button_media.png' alt='Media database' title='Media Database' style="width: 20%;">
  </a>

.. image:: images/buttons/refineGEMs_button_code.png
  :height: 0px
  :width: 0px

.. raw:: html

  <a class="reference external image-reference" href="api.html">
    <img src='_images/refineGEMs_button_code.png' alt='Code' title='Code' style="width: 20%;">
  </a>

Connections to other tools
^^^^^^^^^^^^^^^^^^^^^^^^^^
``refineGEMs`` also provides access points for other tools:

* Reporting the `MEMOTE <https://memote.readthedocs.io/en/latest/index.html>`__ score or generating a complete MEMOTE report
* Updating the SBO-Term annotations based on `SBOannotator <https://github.com/draeger-lab/SBOannotator>`__\ :footcite:p:`Leonidou2023_sboann`
* Balancing the masses and charges with `MassChargeCuration (MCC) <https://github.com/Biomathsys/MassChargeCuration>`__\ :footcite:p:`Mostolizadeh2026_mcc`
* Correcting a biomass objective function with `BOFdat <https://github.com/jclachance/BOFdat>`__

Licence and visual assets
^^^^^^^^^^^^^^^^^^^^^^^^^

The refineGEMs source code is distributed under the MIT licence. The project
logos are licensed separately under `CC BY 4.0
<https://creativecommons.org/licenses/by/4.0/>`__; see the :doc:`brand usage
guide <brand>` for attribution, colours, minimum size, clear space, and
acceptable modifications.

Project-created documentation graphics are also licensed separately under `CC
BY 4.0 <https://creativecommons.org/licenses/by/4.0/>`__; see the
:doc:`documentation graphics licence <graphics>` for covered assets and
third-party exclusions.

Contributing
^^^^^^^^^^^^

Bug reports, feature requests, documentation improvements, and pull requests
are welcome. The repository provides GitHub issue templates; see
:doc:`development` for contribution guidance and review notes.

How to cite
-----------

When using ``refineGEMs``, please cite the latest publication:

  Famke Bäuerle, Gwendolyn O. Döbel, Laura Camus, Simon Heilbronner, and Andreas Dräger.
  Genome-scale metabolic models consistently predict in vitro characteristics of Corynebacterium
  striatum. Front. Bioinform., oct 2023. `doi:10.3389/fbinf.2023.1214074 <https://doi.org/10.3389/fbinf.2023.1214074>`__.



.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   main-modules
   database/intro
   Third-party databases and tools <licenses-citations>
   Content of refineGEMs <api>
   Help and FAQ <help>
   Logo licence and brand usage <brand>
   Documentation graphics licence <graphics>
   Notes for developers <development>

.. footbibliography::
