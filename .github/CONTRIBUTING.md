# Contributing to refineGEMs

## Welcome

Thank you for helping build `refineGEMs`. The project is a Python toolbox for
curating and analysing genome-scale metabolic models, and contributions are
welcome from users, developers, model curators, documentation writers, and
people improving documented modelling workflows.

This guide complements the Sphinx documentation, especially the development,
installation, usage, API, and media database pages on Read the Docs. If this
file and the documentation disagree, please open an issue or include a
documentation fix in your pull request.

## Ways to contribute

Useful contributions include:

- bug reports with reproducible examples;
- fixes for Python code, command line workflows, documentation, examples, and
  packaging;
- new or improved model curation and analysis functionality;
- tests that describe expected behaviour or prevent regressions;
- clearer docstrings, tutorials, and Read the Docs pages;
- corrections to existing media database behaviour or documentation;
- feature ideas, architecture discussions, and reviews of open pull requests.

Small contributions are valuable. If you are unsure whether a change belongs in
the project, open an issue or discussion before investing too much time.

## Reporting Bugs

Please report bugs in the GitHub issue tracker. A good bug report includes:

- the installed `refineGEMs` version or commit;
- the Python version and operating system;
- how `refineGEMs` was installed, including optional dependencies if relevant;
- the command, script, model, medium, or database entry involved;
- the full error message or unexpected output;
- a minimal example that maintainers can run, whenever possible.

Avoid sharing private or unpublished models unless you are allowed to do so. If
the problem depends on data that cannot be shared, describe the relevant model
features, file format, and workflow as precisely as possible.

## Suggesting New Features

Feature requests are easiest to evaluate when they describe the scientific or
workflow problem first. Please include:

- the use case and expected output;
- whether the feature should be part of the Python API, the command line
  interface, the documentation, or several of these;
- examples of related tools, publications, or standards;
- compatibility considerations for existing users and workflows.

For larger changes, start with an issue or discussion so the API, CLI, module
placement, and documentation plan can be agreed before implementation.

## Development Workflow

`refineGEMs` requires Python 3.10 or newer. The latest installation page lists
the currently tested Python versions and optional dependencies.

1. Fork or clone the repository.

   ```console
   git clone https://github.com/draeger-lab/refinegems.git
   cd refinegems
   ```

2. Create and activate a development environment. A conda environment is one
   supported option.

   ```console
   conda create -n refinegems-dev python=3.12
   conda activate refinegems-dev
   ```

3. Install the project in editable mode. Include the documentation extra if you
   will update the Sphinx documentation.

   ```console
   pip install -e ".[docs]" --config-settings editable_mode=strict
   ```

   Some optional workflows need tools such as ModelPolisher, MassChargeCuration,
   BOFdat, or optional extras from `pyproject.toml`. Install only the optional
   dependencies needed for the feature or bug you are working on.

4. Create a feature branch.

   ```console
   git switch -c feature/short-description
   ```

5. Keep changes focused and commit related code, tests, documentation, and
   example updates together.

6. Run the relevant tests and checks before submitting. At minimum, exercise the
   code path you changed. If you add or modify public Python behaviour, add or
   update unit tests. If you modify documentation, build the Sphinx docs.

   ```console
   cd docs
   make html
   ```

## Coding and Documentation Standards

The repository is Python-first. Please keep contributions consistent with the
existing package structure and public interfaces.

- Prefer modular code with clear responsibilities. Add shared helpers only when
  they remove real duplication or clarify an existing workflow.
- Keep the public Python API consistent. Choose names, argument order, return
  values, exceptions, and defaults that match neighbouring functions.
- Preserve backward compatibility where possible. If a breaking change is
  necessary, explain it in the pull request and update documentation and
  examples.
- Follow the existing code style. The repository advertises black formatting;
  format Python changes accordingly.
- Use type hints for new and changed functions where practical.
- Add unit tests for new behaviour, bug fixes, and important edge cases.
- Write docstrings for new public functions, classes, methods, and modules.
  The development documentation uses Google-style docstrings, including the
  local convention of listing arguments as `- name (type): ...` so Sphinx
  renders them cleanly.
- Document new public functionality in the Read the Docs documentation. Public
  Python API changes should appear in API or module pages, and user-facing
  workflows should appear in usage, module, or tutorial pages.
- Keep documentation examples executable. Prefer examples that can be copied
  into a Python session, script, notebook, or shell without hidden setup.
- Add examples when introducing new commands, workflows, or non-obvious
  options.
- Update developer documentation when architecture, module responsibilities, or
  configuration changes.

### Command Line Interface

The command line interface is part of the public user experience.

- Keep command names, option names, argument order, output style, and error
  handling consistent with existing `refinegems` commands.
- Changes to CLI commands require updating both the command help text and the
  Read the Docs documentation.
- Add examples when adding a new command or a new important option.
- Avoid changing existing command behaviour without a compatibility note and a
  migration path where possible.

### Media Database Contributions

A related project is working on separating the media database from the
`refineGEMs` toolbox and building a standalone database with a broader scope,
including connections to different databases and users. Because of that
transition, please coordinate with the maintainers before starting substantial
new media database work in this repository.

Changes that fix current `refineGEMs` media behaviour or documentation should
still be traceable and documented.

- Describe what changed: compound correction, source update, naming change,
  identifier change, documentation-only update, or another scoped maintenance
  change.
- Include the source for the biological or biochemical information, such as a
  publication, database, protocol, or repository file.
- Explain whether the change affects generated media, existing simulations, or
  documented examples.
- Update or regenerate the relevant media and subset documentation pages under
  the Sphinx database documentation.
- When using helper functions such as `Medium.export_to_file(type="docs", ...)`
  or `generate_docs_for_subset(...)`, mention this in the pull request so
  reviewers know how the files were produced.

## Pull Request Checklist

Before opening a pull request, please check that:

- the branch is based on an up-to-date target branch;
- the pull request has a clear title and explains why the change is needed;
- public Python functionality has docstrings and Read the Docs documentation;
- documentation examples are executable;
- new CLI behaviour updates both help text and documentation;
- new commands or workflows include examples;
- relevant unit tests or manual checks were run and are listed in the pull
  request;
- media database maintenance changes were discussed when needed and document
  their sources, scope, and generated documentation updates;
- backward compatibility risks are called out;
- unrelated formatting or generated files are not included.

## Review Process

Maintainers review pull requests for correctness, maintainability, API
consistency, documentation, tests, and compatibility with existing workflows.

Reviews may ask for changes to code structure, names, tests, documentation, or
examples. Please keep discussion focused on the project goals and update the
pull request with follow-up commits. When a review comment is addressed, leave a
short reply if the resolution is not obvious from the diff.

Large or architectural pull requests may need more than one review round. For
major API, CLI, or media-related design changes, discussion before
implementation is strongly encouraged.

## Getting Help

Use GitHub issues for bugs and concrete tasks. Use GitHub discussions or an
issue labelled for discussion when you want feedback on an idea, design,
media-related contribution, or workflow question before opening a pull request.

The Read the Docs pages contain installation instructions, usage examples,
developer notes, API documentation, media database documentation, and help for
known setup issues.

## Citation

If you use `refineGEMs` in scientific work, please cite the latest publication
listed in the repository README and Read the Docs citation page:

Famke Bäuerle, Gwendolyn O. Döbel, Laura Camus, Simon Heilbronner, and Andreas
Dräger. Genome-scale metabolic models consistently predict in vitro
characteristics of Corynebacterium striatum. Frontiers in Bioinformatics, 2023.
https://doi.org/10.3389/fbinf.2023.1214074

Please also cite relevant third-party tools, databases, and publications for
workflows or media database entries that depend on them.
