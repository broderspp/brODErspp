#!/bin/sh

# First of all print some information
echo " "
echo "#### VERIFICATION SCRIPT FOR LARSEN #######"
echo " "
echo "This script will run brODErs++ with the solver SHOCKING."
echo "NOTE THAT THE SOLVER MIGHT NOT STOP! YOU MAY HAVE TO STOP IT"
echo "BY PRESSING CTRL-c, when the simulation is almost over!"
echo " "
echo " "
echo "=========> NOW press RETURN to start.. <========="
echo " "
read

# CTRL-c doesn't stop the script, but only the current running program
trap "echo 'CTRL-C pressed!!!'" 2

# Create outpLL file if it doesn't exist
touch outpSS

# Start reading the file
tail -f outpSS | grep "Sol: " &

# Start Larsen
../../../bin/brODErs++ shocking_input > outpSS

# Grep it.. since there might be a lot of shit in the output
cat outpSS | grep "Sol: " > outpS_tmp

# Kill all the tail processes (still running in background..)
killall tail

# Remove the first 7 characters
cut -c 5- outpS_tmp > shocking_output.dat

# remove the temporary file
rm -f outpS_tmp outpSS

# Run octave
octave octaveResultsPlotter.m
