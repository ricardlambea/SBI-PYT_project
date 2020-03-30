import argparse
from biomabuilderfunctions import *


########################## DEALING WITH ARGUMENTS ##########################
dealer = argparse.ArgumentParser(description =
            "BioMaBuilder is a program for modelling the macro-complex structure of biomolecules formed by proteins and DNA/RNA, starting with each pairing interaction of a complex: protein-protein, protein-DNA/RNA", epilog = " BioMaBuilder. Universidad Pompeu Fabra. MSc Bioinformatics for Health Sciences. 2020. Ricard Lambea-Jane, Carolina Hernandez-Oliver, JV Roig-Genoves")

dealer.add_argument('-i', '--input', required = True, dest = "input_directory", action = "store" , help = "Required argument. It must be a directory with all the PDB files with the pairwise interactions. Review the manual for input file requirements.")

dealer.add_argument ('-fa', '--fasta', required = True, dest = "fasta_file", action = "store", help = "Required argument. It must be a fasta file with all the sequences of the complex")

dealer.add_argument('-o', '--output', required = True, dest = "output_directory", action = "store", help = "Required argument. If the output argument is defined, the program will check out if that directory exists in order to store the PDB file with the final model in it. However, if this directory does not exist, the program will create it and store the PDB file with the final model in it.")

dealer.add_argument('-v', '--verbose', dest = "verbose", action = "store_true", default = False, help= "Optional argument. If applied, the progress of the program will be printed in the standard error.By default is False.")

dealer.add_argument('-sto','--stoichiometry', dest = "stoichiometry", action = "store", default = None, type = int, help = "Optional argument. This argument establishes the number of chains in the final model. It must be a natural number. If not defined, the program will add the number of chains provided in the input directory.")

dealer.add_argument('-rmsd', '--RMSD_threshold', dest = "RMSD_threshold", action = "store", default = 0.6 , type = float, help = "Optional argument. This argument sets up a root-mean-square deviation threshold. By default, the RMSD threshold is 0.6 Angstroms. Review the manual for further information.")

dealer.add_argument('-cs', '--core_selection', dest = "core_selection", action = "store", default = None, help = "Optional argument. This argument allows establishing the reference file that will act as a core for the macro-complex building.")

dealer.add_argument('-ncl', '--number_clashes', dest = "number_clashes", action = "store", default = 30, type = int, help = "Optional argument. This argument sets up the maximum number of clashes allowed during the superimposition. By default, it is defined as 30 clashes. Review the manual for further information.")

dealer.add_argument('-of', '--output_filename', dest = "out_name", action = "store", default = "macrocomplex", help = "Optional argument. This argument sets up the name of the final PDB file. By default, the name of the file will be macrocomplex. It must be defined without any extension (i.e. .pdb)")

args = dealer.parse_args()

if args.input_directory:
    if os.path.isdir(args.input_directory): # Checks out if input_directory argument is a directory.
            pdb_files = obtain_pdb_files(args.input_directory)
            if args.verbose:
                sys.stderr.write("Have been found %d pdb files in %s directory.\n" %(len(pdb_files), args.input_directory))

    else:
        raise IncorrectInput("%s" %(args.input_directory))

############################## CORE OF THE CODE ##############################
if args.core_selection: # Checks out if core_selection argument have been defined.
    core_path =  args.input_directory + "/" + args.core_selection # Uses the file defined in core_selection argument as core file.
    if args.verbose:
        sys.stderr.write("The file %s has been defined as the initial core for macro-complex building.\n" %(args.core_selection))

else: # If core_selection argument has not been defined, obtain the best core file using get_best_core function.
    core_path = args.input_directory + "/" + get_best_core(pdb_files)

    if args.verbose:
        sys.stderr.write("The file %s has been defined as the initial core for macro-complex building.\n" %(core_path))

core_structure = obtain_structure("core", core_path) # Obtains the structure object of the core file.
core_model = core_structure[0] # Obtains the model of the core file.

if args.verbose:
    sys.stderr.write("Starting to construct the macrocomplex...\n")

# Calls the recursive function to construct the final model.
final_model = biobuilder(core_structure = core_structure, files_list = pdb_files, num_iteration = 0, stop_counter = 0, extra_arguments = args)

if args.verbose:
    sys.stderr.write("GREAT! The macrocomplex building process has finished correctly!.\n")

if args.fasta_file: # Aligns fasta sequences
    align_fasta_seqs(args.fasta_file, args.out_name)


if os.path.exists(args.output_directory): # Checks out if output_directory exists
    os.chdir(args.output_directory) # Moves into the output_directory
    if final_model:
        generate_output_file(final_model, args.out_name)
    os.system("mv" + " " + "../biomabuilder/alignments-" + args.out_name + ".txt" + " " + "./" + args.out_name + "-alignments.txt")
    if args.verbose:
        sys.stderr.write("The final model has been stored in %s directory.\n" %(args.output_directory))

else:
    os.mkdir(args.output_directory) # Creates the output_directory
    os.chdir(args.output_directory)
    if final_model:
        generate_output_file(final_model, args.out_name) # Generates the output file with the filename defined.
        os.system("mv" + " " + "../biomabuilder/alignments-" + args.out_name + ".txt" + " " + "./" + args.out_name + "-alignments.txt")
    if args.verbose:
        sys.stderr.write("The final model has been stored in %s directory.\n" %(args.output_directory))
