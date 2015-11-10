#!/bin/bash
#@ wall_clock_limit = 00:20:00
#@ job_name = pos-lulesh-serial
#@ job_type = MPICH
#@ class = test
#@ output = log/pos_lulesh_serial_$(jobid).out
#@ error = log/pos_lulesh_serial_$(jobid).out
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

./bin/lulesh$var > log/$var.log

echo $var >> log/flags.log
echo $(cat log/flags_$var.log) >> log/flags.log
echo $(cat log/$var.log | grep Elapsed) >> log/flags.log
