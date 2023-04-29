# Define the number of atoms and the force constants
num_atoms = 2614
force_constants = [1000, 1000, 1000]

# Open the output file for writing
with open("output.itp", "w") as output_file:

    # Write the header
    output_file.write("; This is a GROMACS topology file\n\n")

    # Write the [ position_restraints ] section header
    output_file.write("[ position_restraints ]\n")

    # Write the column headers
    output_file.write(";  i funct       fcx        fcy        fcz\n")

    # Write the position restraint information for each atom
    for i in range(1, num_atoms+1):
        output_file.write("{0} {1} {2} {3} {4}\n".format(i, 1, *force_constants))

