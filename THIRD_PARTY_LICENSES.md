# Third-party database and tool notices

refineGEMs connects to, queries, and in some cases redistributes selected data
from external resources. This file records notices for redistributed or
licence-sensitive resources. It is not legal advice; users remain responsible
for checking the current upstream terms for their own use case.

The refineGEMs source code is distributed under the MIT licence. Bundled
third-party data, database identifiers, and connected external tools remain
under their own licences or terms.

## BiGG Models

refineGEMs ships BiGG-derived metabolite and reaction tables in its internal
database.

Source: <http://bigg.ucsd.edu/>

Citation:

King ZA, Lu JS, Dräger A, Miller PC, Federowicz S, Lerman JA, Ebrahim A,
Palsson BO, and Lewis NE. BiGG Models: A platform for integrating,
standardizing, and sharing genome-scale models. Nucleic Acids Research
44(D1):D515-D522. <https://doi.org/10.1093/nar/gkv1049>

Licence notice from <http://bigg.ucsd.edu/license>:

> Copyright (c) 2019 The Regents of the University of California
>
> All Rights Reserved
>
> Permission to use, copy, modify and distribute any part of BiGG Models for
> educational, research and non-profit purposes, without fee, and without a
> written agreement is hereby granted, provided that the above copyright notice,
> this paragraph and the following three paragraphs appear in all copies.
>
> Those desiring to incorporate BiGG Models into commercial products or use for
> commercial purposes should contact the Technology Transfer & Intellectual
> Property Services, University of California, San Diego, 9500 Gilman Drive,
> Mail Code 0910, La Jolla, CA 92093-0910, Ph: (858) 534-5815, FAX:
> (858) 534-7345, e-mail: invent@ucsd.edu.
>
> In no event shall the University of California be liable to any party for
> direct, indirect, special, incidental, or consequential damages, including
> lost profits, arising out of the use of this bigg database, even if the
> University of California has been advised of the possibility of such damage.
>
> The BiGG Models provided herein is on an "as is" basis, and the University of
> California has no obligation to provide maintenance, support, updates,
> enhancements, or modifications. The University of California makes no
> representations and extends no warranties of any kind, either implied or
> express, including, but not limited to, the implied warranties of
> merchantability or fitness for a particular purpose, or that the use of the
> BiGG Models will not infringe any patent, trademark or other rights.

## ChEBI

refineGEMs queries ChEBI for compound names, formulae, charges, and identifiers.

Source: <https://www.ebi.ac.uk/chebi/about/>

Citation:

Hastings J, Owen G, Dekker A, Ennis M, Kale N, Muthukrishnan V, Turner S,
Swainston N, Mendes P, Steinbeck C. ChEBI in 2016: Improved services and an
expanding collection of metabolites. Nucleic Acids Research 44(D1):D1214-D1219.
<https://doi.org/10.1093/nar/gkv1031>

ChEBI states that data on its website are available under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) and governed by
EMBL-EBI terms of use.

## KEGG

refineGEMs queries KEGG for genes, enzymes, reactions, compounds, and pathway
annotations. Released refineGEMs files should contain only KEGG identifiers or
annotations, not cached KEGG records or KEGG database exports.

Source: <https://www.kegg.jp/kegg/legal.html>

Citation:

Kanehisa M, Goto S. KEGG: Kyoto Encyclopedia of Genes and Genomes. Nucleic
Acids Research 28(1):27-30. <https://doi.org/10.1093/nar/28.1.27>

KEGG states that academic users may freely use the KEGG website. It also states
that KEGG is not a public database for non-academic use and that non-academic
use requires a commercial licence.

## MetaNetX / MNXref

refineGEMs uses MetaNetX identifiers and mapping tables.

Source: <https://www.metanetx.org/mnxdoc/mnxref.html>

Citation:

Moretti S, Tran VD, Mehl F, Ibberson M, Pagni M. MetaNetX/MNXref: unified
namespace for metabolites and biochemical reactions in the context of metabolic
models. Nucleic Acids Research 49(D1):D570-D574.
<https://doi.org/10.1093/nar/gkaa992>

MetaNetX states that, except where otherwise noted, data from the site are
licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
MetaNetX also notes that MNXref data can be subject to the original licence
restrictions of external source resources.

## ModelSEED

refineGEMs ships ModelSEED-derived compound data in its internal database.

Source: <https://github.com/ModelSEED/ModelSEEDDatabase>

Citation:

Henry CS, DeJongh M, Best AA, Frybarger PM, Linsay B, Stevens RL.
High-throughput generation, optimization and analysis of genome-scale metabolic
models. Nature Biotechnology 28:977-982.
<https://doi.org/10.1038/nbt.1672>

The ModelSEEDDatabase licence states that, except where otherwise noted, data in
the ModelSEED GitHub repository are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). It also states that
ModelSEED includes information sourced in part from external resources, and that
data derived from those resources may remain subject to the original resource
restrictions.

## Adapted MEMOTE code

Parts of `src/refinegems/curation/biomass.py` were adapted from MEMOTE and
modified by Gwendolyn O. Döbel. MEMOTE is distributed under the Apache License
2.0.

Source: <https://github.com/opencobra/memote>

Citation:

Lieven C, Beber ME, Olivier BG, et al. MEMOTE for standardized genome-scale
metabolic model testing. Nature Biotechnology 38:272-276.
<https://doi.org/10.1038/s41587-020-0446-y>

Licence: <https://github.com/opencobra/memote/blob/develop/LICENSE>

## Connected tools

The following connected software tools are distributed by their upstream
projects under their own software licences:

- BOFdat: Original software under [MIT](https://choosealicense.com/licenses/mit/), used fork by Draeger lab, also under
  MIT; cite <https://doi.org/10.1371/journal.pcbi.1006971>.
- MassChargeCuration (MCC): [GNU LGPL-3.0](https://choosealicense.com/licenses/lgpl-3.0/); cite <https://doi.org/10.1128/spectrum.03200-24>.
- MEMOTE: [Apache-2.0](https://choosealicense.com/licenses/apache-2.0/); cite <https://doi.org/10.1038/s41587-020-0446-y>.
- SBOannotator: [GNU GPL-3.0](https://choosealicense.com/licenses/gpl-3.0/); cite <https://doi.org/10.1093/bioinformatics/btad437>.
