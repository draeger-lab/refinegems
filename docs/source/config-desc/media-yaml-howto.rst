How to adjust the media config file
===================================

The media configuration file is divided into two main parts:

- params: Sets the defaults for all non-set parameters in the media definitions
- media: Actual part where the media are listed

Media config: params block
--------------------------

This block has a total of four options that can be set. If not given, uses a default:

Example:

.. code:: yaml

    params:
        aerobic: True 
        supplement: None
        default_flux: 10.0
        o2_percent: 1.0

Options:

- ``aerobic`` boolean

    Sets the media to aerobic, if True (adds O2 when necessary).
    Default sets the media to anaerobic.

- ``supplement`` None, std or min

    Set the supplement mode. If std/min, uses the minimal medium/standard exchanges to 
    supplement the medium to enable growth, if growth is not possible.

- ``default_flux`` float

    Default flux to set if no flux value is given for a substance. Defaults to 10.0.

- ``o2_percent`` float

    Percentage of flux to use for oxygen. Defaults to 1.0, which is 100%, meaning the normal (default) flux.



Media config: media block
-------------------------

The media block lists the actual media plus their modifications to be loaded.
Its structure is as follows:

.. code:: yaml

    media:
        name1:
            option1: value1
            option2:
                value21: mappedvalue1
                value22: mappedvalue2
        name 2:
        ...

The block begins with the keyword ``media``. This block can contain an arbitrary number of subdictionaries. Each starts 
with a name. The placeholder ``name`` can be substituted with a name for a medium or a name of a medium in the media 
database. If the ``name`` of the medium is not a valid name for the media database, either base or external_base is 
required to set the base medium for the medium ``name``. The following options can be used to further modify each medium:

.. note::

    To specify the medium ``name`` without modification, simply set the name and leave the dictionary empty.

- ``base`` None or str

     The value for base has to be a valid name/abbreviation of the media database.

- ``external_base`` None or str(path)

    Similar to base, but expects a path to an external medium file.

- ``add_subset`` List of subset name: percentage (float, e.g. 0.5)

    Add subsets of the media database to the medium. Requires the abbreviation of the 
    subset and a percentage of how much in relation to the base should be added.

- ``add_medium`` List of medium name: mode of combination

    Similar to add_subset but adds another medium instead of a subset. The mode of
    combination can be:

    - a single float, e.g. ``0.1``; 
    - a tuple of two floats, e.g. ``(0.1, 0.2)``; 
    - ``"+"`` to simply add the media together; 
    - or ``NULL`` to set the fluxes to ``None``.

- ``add_external`` List of paths: percentage (float, e.g. 0.5)

    Similar to add_subset, but adds a medium loaded from an external path and
    adds it to the medium.

- ``add_substance`` Name:  Null|float|str in format 'X.X%' with X being [0-9]*

    Add more substances from the database to the medium or change the flux values of 
    already existing ones. This has the highest priority, overwriting all other flux values.
    When a percentage is given uses the most specific set default (medium > params) as the reference point.

.. hint::

    In addition to these options, all options under params can be set here separately
    for each medium as well, overwriting the default only for that specific medium.

Example for aerobic LB with lower fluxes, but double oxygen and added glycerol:


.. code:: yaml

    ...
        LB_mod:
            aerobic: True
            default_flux: 5.0
            o2_percent: 2.0
            add_substance:
                Glycerol: 0.87
    ...
