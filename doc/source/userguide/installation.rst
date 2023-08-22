=============
Installation
=============

Using anaconda
--------------------------


Anaconda
"""""""""

Creating environments using Anaconda is recommended due to the ability to
create more than one environment. It is also recommended because you can
keep dependencies separate from one another that might conflict if you had
them all in your root environment. For example, if you had all the dependencies
for a Pandas environment and all the dependencies for a Pyrad environment in
your root environment, there might be conflicts between channels and packages.
So Anaconda allows you to create multiple environments to avoid these issues.

To download and install `Anaconda <https://www.anaconda.com/download/#>`_.

While Anaconda is downloading, it will ask if you want to set a path to it, or
let Anaconda set a default path. After choosing, Anaconda should finish
downloading. After it is done, exit the terminal and open a new one to make
sure the environment path is set. If conda command is not found, there is help
on running conda and fixing the environment path, found here:

* `How to Run Conda <https://stackoverflow.com/questions/18675907/how-to-run-conda>`_

Setting a Channel
"""""""""""""""""""

Anaconda has a cloud that stores many of its packages. It is recommended, at
times, to use the conda-forge channel instead. Conda-Forge is a community led
collection of packages, and typically contains the most recent versions of the
packages required for Pyrad. Also Pyrad is on Conda-Forge. Having packages in
an environment, within the same channel, helps avoid conflict issues. To add
conda-forge as the priority channel, simply do::

        conda config --add channels conda-forge

You can also just flag the channel when conda install packages such as::

        conda install -c conda-forge numpy

More on managing channels can be found here:

* `Managing Channels <https://conda.io/docs/user-guide/tasks/manage-channels.html>`_

.. _condaenv:

Creating an Environment
"""""""""""""""""""""""""

There are a few ways to create a conda environment for using Pyrad or other
packages. One way is to use the environment file, found here:

* https://github.com/MeteoSwiss/pyrad/blob/master/environment.yml

To create an environment using this file, use the command::

        conda env create -f environment.yml

This will then create an environment called pyrad_env that can be activated
by::

        conda activate pyrad_env

or deactivated after use::

        conda deactivate pyrad_env

Once the environment is created and activated, you can install more packages
into the environment by simply conda installing them. An example of this is,
if you want Jupyter Notebook to run in that enviroment with those packages::

        conda install -c conda-forge jupyter notebook

while that environment is activated. Another way to create a conda environment
is by doing it from scratch using the conda create command. An example of this::

        conda create -n pyrad_env -c conda-forge python=3.10 pyrad_mch pyart_mch netCDF4
        cartopy scipy numpy matplotlib

This will also create an environment called pyrad_env that can be activate the
same way, as mentioned above. To then run your coding editor within the
environment, run in the command line::

        python

or::

        ipython

or::

        jupyter notebook

or even::

        spyder

depending on what you installed in your environment and want to use for coding.

Environment variables
"""""""""""""""""""""""""

Pyrad uses some environment variables to access configuration and shared library files. These environment variables are the following:

PYART_CONFIG
    Path to the Py-ART configuration file, for example at MeteoSwiss, we use the `following one <https://github.com/MeteoSwiss/pyrad/blob/master/config/pyart/mch_config.py>`_.

METRANETLIB_PATH
    Path to the metranet shared library files, which are required to read the Metranet format of the Swiss operational radars.
    If you are interest in getting access to these files, please contact daniel.wolfensberger@meteoswiss.ch.

RSL_PATH
    Path to the RSL library which is used by some routines of Py-ART. Please be aware that this library is deprecated, and the Py-ART is currently
    working on an alternative.

A good way to avoid redefining these environment variables everytime you want to use Pyrad/Py-ART, is to define them in your .bashrc file.


More Information
"""""""""""""""""""""""""


For more an conda and help with conda:

* https://conda.io/docs/
* https://gitter.im/conda/conda


Using pip
--------------------------

Another possibility is to use the pip package manager to install pyrad and its dependencies::

        python3 -m pip install pyart_mch pyrad_mch numpy scipy matplotlib netcdf4 xarray trmm_rs cartopy


Another way is to use the requirements.txt file, found here:

* https://github.com/MeteoSwiss/pyrad/blob/master/requirements.txt

To create an environment using this file, use the command::

        python3 -m pip install -r requirements.txt

Before installing Pyrad with pip it is highly recommended to first create a virtual environment either with conda :ref:`condaenv` or with `pip <https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/>`_

From the source
--------------------------

Getting the code
"""""""""""""""""""""""""

To get a copy of the Pyrad superproject simply place yourself in the  desired working directory (It is strongly recommended to use your $HOME in order to be able to use 
some of the Pyrad tools) and type::

        git clone --recursive https://github.com/MeteoSwiss/pyrad.git 

The recursive keyword fetches automatically all the submodules depending on the main superproject.

Regular users should use the “master” branches of both Pyrad and Py-ART. To check that you use 
the “master” branch of Pyrad place yourself in the root directory of the project and type::
        git branch

And eventually::

        git checkout master

And to check that you use the “master” branch of Py-ART go to the directory src/pyart and repeat the 
procedure above

MeteoSwiss developers should use instead the “dev” branch for both Pyrad and Py-ART. PyTDA only 
has a master branch.

Compilation
"""""""""""""""""""""""""

For the initial compilation of the software activate the conda environment, i.e.::
        conda activate pyrad

Then go to pyrad/src and execute::
        make_all.sh

This command takes care of compiling, Py-ART, PyTDA and Pyrad. To compile them separately you 
can use the scripts make_pyart.sh, make_pytda.sh and make_pyrad.sh or compile them separately by moving to all subdirectory of pyrad/src and run::

        python -m pip install . 


Py-ART has a default config file called default_config.py located in folder pyart. If you would like to 
work with a different config file you have to specify the location in the variable PYART_CONFIG in 
your conda environment file. 

For example::
        export PYART_CONFIG= [Pyrad_path]/config/pyart/mch_config.py

The Pyrad library has its own config file in the aforementioned path


