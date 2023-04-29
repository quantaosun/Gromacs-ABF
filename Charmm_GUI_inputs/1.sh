#!/bin/bash

# Define the number of atoms and the force constants
num_atoms=2614
force_constants="1000 1000 1000"

# Open the output file for writing
echo "[ position_restraints ]" > output.itp
echo ";  i funct       fcx        fcy        fcz" >> output.itp

# Write the position restraint information for each atom
for i in $(seq 1 $num_atoms)
do
    echo "$i   1   $force_constants" >> output.itp
done

