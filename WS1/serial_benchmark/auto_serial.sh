#!/bin/bash
rm -rf ../lulesh/log/*
rm -rf ../lulesh/bin/*
int2array ()
{
  input=$1;
  local -a array;
  i=0;

  while [ $input != 0 ]
  do
    array[$i]=$[$input%10]
    input=$[$input/10];

    i=$[$i + 1];
  done
  return=${array[@]};
}

flags_len=7; # Number of flags.
declare -a flags=("-march=native" "-fomit-frame-pointer" "-floop-block" "-floop-interchange" "-floop-strip-mine" "-funroll-loops" "-flto") # Array with flags.

iteration_cnt=1; # Counter of iterations in the most outer loop.

for i in $(seq 1 $((2**($flags_len)))); # Let's iterate over all of possible combinations.
do
  # Convert int(i) to binary(i)
  binary_number=$(echo "obase=2;$i" | bc);

  # Convert binary to array of integers
  int2array $binary_number
  array=${return[@]};

  status_cnt=0; # counter for the inner loop.
  declare -a current_flags=(); # string containing all of the flags that have to ba added to Makefile.

  # Check which flags are turned on in the current combination and add them to the current_flags array.
  for status in ${array[@]};
  do
    if [ $status -eq 1 ]
      then current_flags="$current_flags ${flags[$status_cnt]}"
    fi
    status_cnt=$[$status_cnt + 1]
  done

  # Until here, the i-th combination is set up.

  # Let's make our Makefile default.
  cp ../lulesh/Makefile_serial ../lulesh/Makefile
  # Let's write the new flag combination to Makefile.
  sed -i "22i\\CXXFLAGS= -O3${current_flags}" ../lulesh/Makefile
  sed -i "23i\\LDFLAGS= -O3${current_flags}" ../lulesh/Makefile
  # Make binary with unique id belonging to current combination.
  sed -i "6i\\LULESH_EXEC = bin/lulesh${iteration_cnt}" ../lulesh/Makefile

  # Write the unique ID directly to ll_serial.sh that you're posting.
  cp ../serial_benchmark/ll_serial_flags.sh ../serial_benchmark/ll_serial_flags_${iteration_cnt}.sh
  sed -i "19i\\var=${iteration_cnt}" ../serial_benchmark/ll_serial_flags_${iteration_cnt}.sh

  # Rebuild the project.
  cd ../lulesh
  make clean && make

  ### RUN AT HOME (comment out if you run it on SUPERMUC!)
  # Create a log file and run the lulesh.
  ./bin/lulesh${iteration_cnt} > log/$iteration_cnt.log
  # Which flags (with id) belong to which measure time?
  echo $iteration_cnt >> log/flags.log
  echo $current_flags >> log/flags.log
  echo $(cat log/$iteration_cnt.log | grep Elapsed) >> log/flags.log

  ### RUN ON SUPERMUC (comment out if you run it at HOME!)
  #echo $current_flags >> log/flags_$iteration_cnt.log
  #llsubmit ../serial_benchmark/ll_serial_flags_${iteration_cnt}.sh

  iteration_cnt=$[iteration_cnt + 1];
  echo $iteration_cnt
done
