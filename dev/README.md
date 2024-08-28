# Developer Zone

This is a directory containing information, notebooks with open work and testing examples and much more for the developers of the `refineGEMs` package. It also serves as an additional hub, where developers can share and discuss their code.

Additionally, some helpful code for making developement easier and standardised 
between developers has been added here as well.

More information about the different files collected here can be found below.

## Help with Coding

- `docstring-format.mustache`<br>
  As the name suggest, this file contains the basic setup for the docstring
  format used for this package. It can be integrated into your IDE (e.g. VSCode) to allow autodoc to use the correct style.
- `TodoTree_params.txt`<br>
  If you use VSCode, you can install the TodoTree extension and add this file to your settings.json. It will highlight all listed keyword can collect them in a seperate tab. This makes it easier to track TODOs, points of DISCUSSIONs and more and allows faster identification of open issues inside the code.

## Working Progress Notebooks

- `db_extension` : <br>
  Functions for adding new information to the media database. Should be added to either database or medium modules in the future.<br>
  Status: Almost complete
- `under_construction` :<br>
  Place to collect code snippets, ideas and more that have no playce (or not yet) anywhere else in the codespace.<br>
  Status: Hell, but feel free to visit<br>
- `growth_curves` : <br>
  Ideas for reading in growth plate data and calculating growth curves from it.<br>
  Status: Started, but stuck
- `gapfill_testing`:<br>
  Blocks for testing the different gap fill algorithmns<br>
  Status: Debugging  

## Ideas

- `idea_gapfill.png`: <br>
  Idea of how to revamp the gap fill model. <br>
  Statius: Implemented

## Removed from Package

- `old_gapfill.rst`: <br>
  Old documentation of the refineGEMs pre-2.0.0 gapfill module. <br>
  Status: Re-check for relevance and delete.
