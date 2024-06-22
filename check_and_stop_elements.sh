#!/bin/bash

ArrayID=$1

while true; do
  numOfRunning=$(squeue -l -j $ArrayID | grep RUNNING | wc -l)
  echo
  echo "Pending elements in array job $ArrayID:"
  squeue -l -j $ArrayID | grep PENDING
  echo
  echo "Number of running elements in array job $ArrayID:"
  echo "$numOfRunning"
  echo
  echo "Checking the stop string \"Running model model_1_multimer_v3_pred_0\" in alphaMSAs.*.err for each running element:"
  echo

  for j in $(squeue -l -j $ArrayID | grep RUNNING | awk '{print $1}')
  do
    jobid=$(scontrol show job $j | grep JobId= | awk '{print $1}' | sed -nr 's|JobId=()|\1|p')
    echo -n "Checking alphaMSAs.${jobid}.err ... "
    grep -m 1 "Running model model_1_multimer_v3_pred_0" alphaMSAs.${jobid}.err 2>&1 1>/dev/null
    if [ $? -eq 0 ]; then
      echo -en "\e[31mthe stop string is found!\e[0m ... go to cancel job ${jobid} ($j) ... scancel $j ... \e[32mjob canceled.\e[0m\n"
      scancel $j
    else
      echo -e "the stop string is NOT found!"
    fi
   done
   sleep 300
   clear
done
