#!/bin/bash
#@ wall_clock_limit = 00:20:00
#@ job_name = pos-lulesh-serial
#@ job_type = serial
#@ class = test
#@ output = log/pos_lulesh_hybrid_$(jobid).out
#@ error = log/pos_lulesh_hybrid_$(jobid).out
#@ node = 1
#@ total_tasks = 1
#@ node_usage = not_shared
#@ energy_policy_tag = lulesh
#@ minimize_time_to_solution = yes
#@ island_count = 1
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh
. $HOME/.bashrc

./lulesh2.0
