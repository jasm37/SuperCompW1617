#!/bin/bash
#@ wall_clock_limit = 00:20:00
#@ job_name = pos-amg2013-mpi
#@ job_type = MPICH
#@ class = test
#@ output = pos_amg2013_mpi_16_$(jobid).out
#@ error = pos_amg2013_mpi_16_$(jobid).out
#@ node = 1
#@ total_tasks = 16
#@ node_usage = not_shared
#@ energy_policy_tag = amg2013
#@ minimize_time_to_solution = yes
#@ island_count = 1
#@ notification = never
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module unload mpi.ibm
module load mpi.intel

mpiexec -n 16 ./amg2013  -laplace -n 160 160 160 -P 4 2 2 -printstats
