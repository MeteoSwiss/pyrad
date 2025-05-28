=============
Debugging
=============

Debugging of pyrad can be done by setting the environment variable **PYRAD_DEBUG** to 1::

   export PYRAD_DEBUG=1

Then when calling *main_process_data.py* or any other pyrad scripts, the *pdb* debugger will start automatically if a warning is encountered. You can then freely interact with the current state of variables at the moment of the warning.
