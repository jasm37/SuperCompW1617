#!/bin/bash
#@ wall_clock_limit = 00:20:00
#@ job_name = pos-amg2013-mpi
#@ job_type = MPICH
#@ class = test
#@ output = pos_amg2013_mpi_2_unroll_simd.out
#@ error = pos_amg2013_mpi_2_unroll.out
#@ node = 1
#@ total_tasks = 2
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

mpiexec -n 2 ./amg2013  -laplace -n 320 320 320 -P 2 1 1 -printstats
