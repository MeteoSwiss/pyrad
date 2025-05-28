========================
Contributing to Pyrad
========================

Pyrad is an open project and all contributions are welcome.

Please check our `User manual <https://github.com/MeteoSwiss/pyrad/blob/master/additional_doc/pyrad_user_manual.pdf>`_ to get started!

Code quality
"""""""""""""""""""""""""

The code of pyrad is quality-checked with `ruff <https://docs.astral.sh/ruff/>`_ (PEP8) and `black <https://github.com/psf/black>`_ (code formatting).
Any code that does not follow these guidelines will make the CI of the commit or PR fail. The best way to avoid this is to install pre-commit hooks that will run ruff and black for every local commit you make.

To install the pre-commit hooks simply go to the root directory or your pyrad repository (the one with the .git directory) and run::

    pip install pre-commit
    pre-commit install
