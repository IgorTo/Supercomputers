#!/bin/bash
#@ wall_clock_limit = 00:20:00
#@ job_name = pos-lulesh-hybrid
#@ job_type = MPICH
#@ class = test
#@ output = log/hybrid_8_2.out
#@ error = log/pos_lulesh_hybrid_$(jobid).out
#@ node = 1
#@ total_tasks = 16
#@ node_usage = not_shared
#@ energy_policy_tag = lulesh
#@ minimize_time_to_solution = yes
#@ island_count = 1
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh
. $HOME/.bashrc

module unload mpi.ibm
module load mpi.intel

export OMP_NUM_THREADS=2 # Correct here for different # of threads.
mpiexec -n 8 ./lulesh2.0
