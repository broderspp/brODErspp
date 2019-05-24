#!/bin/sh

# First of all print some information
echo " "
echo "#### VERIFICATION SCRIPT FOR LARSEN #######"
echo " "
echo "This script will run brODErs++ with the solver LARSEN."
echo "NOTE THAT THE SOLVER MIGHT NOT STOP! YOU MAY HAVE TO STOP "
echo "IT BY PRESSING CTRL-c, when the simulation is almost over!"
echo " "
echo " "
echo "=========> NOW press RETURN to start.. <========="
read

# CTRL-c doesn't stop the script, but only the current running program
trap "echo 'CTRL-C pressed!!!'" 2

# Create outpLL file if it doesn't exist
touch outpLL

# Start reading the file
tail -f outpLL | grep "Sol: " &

# Start Larsen
../../../bin/brODErs++ larsen_input.in > outpLL

# Grep it.. since there might be a lot of shit in the output
cat outpLL | grep "Sol: " > outpL_tmp

# Kill all the tail processes (still running in background..)
killall tail

# Remove the first 7 characters
cut -c 5- outpL_tmp > outpL

# remove the temporary file
rm -f outpL_tmp outpLL

# Run octave
octave octaveResultsPlotter.m
