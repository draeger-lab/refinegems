Third-party databases and tools
===============================

``refineGEMs`` connects to, queries, and in some cases redistributes selected
data from external databases and tools. This page summarises the expected
citations and the license or terms that are relevant for use of these resources.
Redistribution notices are also collected in the repository-level
``THIRD_PARTY_LICENSES.md`` file.

.. important::

    This page is a documentation aid for users and maintainers, not legal
    advice. Users remain responsible for checking the current upstream terms,
    especially for commercial use or redistribution of downloaded database
    content.

:icon:`fa-solid fa-database` Databases
--------------------------------------

.. list-table::
   :class: no-scrollbar-table
   :header-rows: 1
   :widths: 20 25 30 25

   * - Resource
     - Use in ``refineGEMs``
     - License or terms
     - Citation
   * - BiGG Models
     - Local BiGG metabolite and reaction tables are shipped in the internal
       database and BiGG identifiers are used for media export, mapping, and
       model curation.
     - Free for educational, research, and non-profit purposes if the BiGG
       copyright notice and license text are retained. Commercial use requires
       a separate agreement with UC San Diego.
     - :footcite:t:`King2016_bigg`
   * - ChEBI
     - ChEBI is queried for compound names, formulae, charges, and identifiers.
     - ChEBI data are available under `CC BY 4.0
       <https://creativecommons.org/licenses/by/4.0/>`__ and governed by
       EMBL-EBI terms of use.
     - :footcite:t:`Hastings2016_chebi`
   * - KEGG
     - KEGG is queried for genes, enzymes, reactions, compounds, and pathway
       annotations. Released ``refineGEMs`` files should contain only KEGG
       identifiers or annotations, not cached KEGG records or KEGG database
       exports.
     - Academic users may use the KEGG website freely. KEGG is not a public
       database for non-academic use; non-academic use and service provision
       require the appropriate KEGG license.
     - :footcite:t:`Kanehisa2000_kegg`
   * - MetaNetX / MNXref
     - MetaNetX identifiers and mapping tables are used for reaction and
       metabolite mapping, including gap-filling.
     - MetaNetX data are `CC BY 4.0
       <https://creativecommons.org/licenses/by/4.0/>`__ unless otherwise
       noted. MNXref records sourced from external databases may remain subject
       to the original database restrictions.
     - :footcite:t:`Moretti2021_metanetx`
   * - ModelSEED
     - Local ModelSEED compound data are shipped in the internal database and
       used for charge and formula comparison.
     - ModelSEED data are `CC BY 4.0
       <https://creativecommons.org/licenses/by/4.0/>`__ unless otherwise
       noted. Records derived from external resources may remain subject to the
       original resource restrictions.
     - :footcite:t:`Henry2010_modelseed`

:icon:`fa-solid fa-screwdriver-wrench` Tools
--------------------------------------------

.. list-table::
   :class: no-scrollbar-table
   :header-rows: 1
   :widths: 20 25 30 25

   * - Tool
     - Use in ``refineGEMs``
     - License or terms
     - Citation
   * - BOFdat
     - Called to generate selected biomass objective function coefficients from
       experimental input data.
     - Distributed as open-source software under the `MIT license <https://choosealicense.com/licenses/mit/>`__; check the installed upstream or
       forked package for the exact license text.
     - :footcite:t:`Lachance2019_bofdat`
   * - MassChargeCuration (MCC)
     - Called to curate metabolite masses and charges.
     - Distributed as free software under the `GNU LGPL-3.0 license <https://choosealicense.com/licenses/lgpl-3.0/>`__.
     - :footcite:t:`Mostolizadeh2026_mcc`
   * - MEMOTE
     - Used to run model consistency tests and generate MEMOTE reports. Some
       consistency helper functions are imported for internal checks, and
       biomass-normalisation code in ``refineGEMs`` was adapted from MEMOTE.
     - Distributed under the `Apache-2.0 license <https://choosealicense.com/licenses/apache-2.0/>`__.
     - :footcite:t:`Lieven2020_memote`
   * - SBOannotator
     - Called to assign SBO terms to SBML model entities.
     - Distributed as free software under the `GNU GPL-3.0 license
       <https://choosealicense.com/licenses/gpl-3.0/>`__. This affects
       users who install or redistribute the optional SBOannotator-dependent
       functionality.
     - :footcite:t:`Leonidou2023_sboann`

:icon:`fa-solid fa-file-shield` Compliance notes for maintainers
-----------------------------------------------------------------

The internal SQLite database currently redistributes BiGG and ModelSEED-derived
tables. When these tables are updated, preserve the upstream license notices,
record the source version or download date, and keep this page in sync.

The ``refineGEMs`` source code is distributed under the MIT license. Bundled
third-party data, database identifiers, and connected external tools remain
under their own licenses or terms, as listed above and in
``THIRD_PARTY_LICENSES.md``.

KEGG-dependent functionality should be documented and implemented as online
queries through the user-accessed KEGG services. Avoid bundling KEGG data in
released artifacts unless the project has confirmed that the intended release
and user base are covered by an appropriate KEGG license.

Because ModelSEED and MetaNetX aggregate data from other resources, downstream
use may need to follow the original source terms for specific records. Preserve
provenance columns and identifiers where possible so users can trace records
back to their source database.

.. footbibliography::
