#!/bin/bash
#@ wall_clock_limit = 00:20:00
#@ job_name = pos-amg2013-hybrid
#@ job_type = MPICH
#@ class = test
#@ output = pos_amg2013_hybrid_4x4_$(jobid).out
#@ error = pos_amg2013_hybrid_4x4_$(jobid).out
#@ node = 1
#@ total_tasks = 4
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

export OMP_NUM_THREADS=4
mpiexec -n 4 ./amg2013  -laplace -n 320 320 160 -P 2 1 2 -printstats
