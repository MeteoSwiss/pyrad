#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
================================================
generate_cartesian_tasks_file
================================================

This script generates a tasks file for parallel processing of rad4alp
Cartesian data at the CSCS.
DO NOT MODIFY!!! To use it make a local copy

"""

# Author: fvj
# License: BSD 3 clause

import datetime
import atexit
import os
import pandas as pd
import numpy as np

print(__doc__)


def main():
    """
    """
    print("====== tasks file generation started: %s" %
          datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))
    atexit.register(_print_end_msg,
                    "====== tasks file generation finished: ")
    fpath = os.path.expanduser('~') + '/pyrad/tools/processData/'
    tasks_name = 'rad4alp_gc_freq'

    ftasks_name = fpath + 'tasks_' + tasks_name
    task = (
        os.path.expanduser('~') +
        '/pyrad/tools/processData/get_and_process_rad4alp_cartesian_data_cscs.sh')

    # generate list of dates
    start_date = datetime.datetime(2020, 7, 23)
    ndays = 1
    datelist = pd.date_range(start_date, periods=ndays).tolist()

    # datelist.append(datetime.datetime(2018, 3, 7))

    # create launcher
    # KESCH specs:
    # nodes: 12
    # CPUs per node: 2
    # Cores per CPU: 12
    # Total cores per node : 24
    # Memory per node 256 GB
    # Total processors : 288
    # Total memory : 3 TB
    #
    # To check specs of partition use: scontrol show partition
    # Partition normal:
    #   nodes : 11
    #   CPUs (processors): 264
    #   Maximum processing time: 1 day
    #
    # Partition postproc:
    #   nodes : 5
    #   CPUs (processors) : 120
    #   Maximum processing time : 1 day
    ntasks = ndays
    partition = 'postproc'
    cpu_time = '24:00:00'

    # Parallelization commands. Set to 1 if you want it true
    MP_DSET = 0
    MP_PROD = 0
    MP_PROF = 0

    # add CPUs per task if task is parallelized
    cpus_per_task = 1

    cpus_per_node = 24
    if partition == 'postproc':
        max_cpus = 120
    elif partition == 'normal':
        max_cpus = 264

    cpus_t = int(ntasks * cpus_per_task)
    if cpus_t > max_cpus:
        warn('Number of CPUs to use to large. Available=' + str(max_cpus) +
             ' Required=' + str(cpus_t))
        cpus_t = max_cpus

    print('ntasks: ' + str(ntasks))
    print('cpus_t: ' + str(cpus_t))

    nodes = int(np.ceil(cpus_t / cpus_per_node))
    ntasks_per_node = int(np.ceil(ntasks / nodes))
    ntasks_per_core = int(np.ceil(ntasks_per_node / cpus_per_node))

    flauncher_name = fpath + 'launch_greasy_' + tasks_name + '.sbatch'
    with open(flauncher_name, 'w', newline='') as txtfile:
        txtfile.write(
            '#!/bin/bash -l\n' +
            '#SBATCH --time=' + cpu_time + '\n' +
            '#SBATCH --nodes=' + str(nodes) + '\n' +
            '#SBATCH --ntasks-per-node=' + str(ntasks_per_node) + '\n' +
            '#SBATCH --ntasks-per-core=' + str(ntasks_per_core) + '\n' +
            '#SBATCH --cpus-per-task=' + str(cpus_per_task) + '\n' +
            '#SBATCH --partition=' + str(partition) + '\n' +
            '#SBATCH --account=msrad\n\n' +
            'module use /users/jfigui/easybuild/kesch/modules/all\n' +
            'module load greasy/2.1-cscs-gmvolf-17.02\n\n' +
            'export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK\n' +
            'greasy -f ' + ftasks_name)

        print('written launcher ' + flauncher_name)

    # start and end time of processing each date
    start_time = '000001'
    end_time = '240000'

    # if more than one parameter separate by coma
    # info has to have the same number of cfgfiles to be applied
    # cfgpath can have one value or as many as cfg files
    cfgfiles = (
        'cscs_rad4alp_POH.txt')
    info = 'POH'
    cfgpath = os.path.expanduser('~') + '/pyrad/config/processing/'

    # Separate products to be retrieved with commas
    get_data = 1
    rm_data = 1
    prod = 'RZC'
    data_destbase = '/store_new/mch/msrad/radar/rad4alp/rawdata/'

    with open(ftasks_name, 'w', newline='') as txtfile:
        for day in datelist:
            txtfile.write(
                task + ' -d ' + day.strftime("%Y%m%d") +
                ' --start_time ' + start_time + ' --end_time ' + end_time +
                ' -c ' + cfgfiles + ' --prod ' + prod + ' -i ' + info +
                ' --cfgpath ' + cfgpath +
                ' --get_data ' + str(get_data) + ' --rm_data ' + str(rm_data) +
                ' --data_destbase ' + data_destbase +
                ' --MP_DSET ' + str(MP_DSET) + ' --MP_PROD ' + str(MP_PROD) +
                ' --MP_PROF ' + str(MP_PROF) + '\n')

        txtfile.close()

        print('written task file ' + ftasks_name)


def _print_end_msg(text):
    """
    prints end message

    Parameters
    ----------
    text : str
        the text to be printed

    Returns
    -------
    Nothing

    """
    print(text + datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M:%S"))


# ---------------------------------------------------------
# Start main:
# ---------------------------------------------------------
if __name__ == "__main__":
    main()
