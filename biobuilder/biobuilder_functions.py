import sys
import os
import re
from Bio.PDB import PDBParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB import Superimposer
from Bio.PDB import NeighborSearch

class IncorrectInput(NameError):
	def __init__(self,input):
		self.input = input
	def __str__(self):
		return "The directory name %s is not a directory " %(self.input)

def obtain_pdb_files (directory):
	""" This function obtains the ".pdb" files from a specified directory and returns a list of the aforementioned files. """

	directory = os.listdir(directory) # Saves all the files in a variable.
	list_of_pdb_files = []
	for file in directory:
		if file.endswith(".pdb"):
			list_of_pdb_files.append(file) # Saves the pdb files in a list.
	return list_of_pdb_files


def obtain_structure(pdb_name,pdb_file):
    """ This function parses a PDB file to obtain a structure object. """

    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdb_file[0:-4], pdb_file)
    return structure

def get_backbone_atoms_protein (chain):
	""" This function returns a list of proteins' backbone atoms (CA) from a given chain. """

	ca_atoms = []
	for residue in chain:
		if residue.get_id()[0] == " " and residue.has_id("CA"):
			ca_atoms.append(residue['CA'])
	return ca_atoms

def get_backbone_atoms_nucleicacids (chain):
	""" This function returns a list of nucleic acids' backbone atoms (C4') from a given chain. """

	c4_atoms = []
	for residue in chain:
		if residue.get_id()[0] == " " and residue.has_id("C4\'"):
			c4_atoms.append(residue["C4\'"])
	return c4_atoms

def get_molecule_type (chain):
	""" This function identifies the molecule type (RNA, DNA, PROTEIN) of a given chain and returns it as a string. """
	# Creates a list for each DNA and RNA chain with all possible letter for eachone.
	RNA = ['A','U','C','G','I']
	DNA = ['DA','DT','DC','DG','DI']
	molecule_type = ""

	for residue in chain:
		residue_name = residue.get_resname().strip()
		break
	if residue_name in RNA:
		molecule_type = "RNA"
	elif residue_name in DNA:
		molecule_type = "DNA"
	else:
		molecule_type = "Protein"
	return molecule_type


def get_best_core(list_of_pdb_files):
	""" This function selects the file containing the chain with larger number of interactions from the list of input files to act as the core for the superimposition strategy, and returns this file back. """

	recount = {}
	## Creates a dictionary with the different chains in input_directory files and count the times they appear ##
	for file in list_of_pdb_files:
		if file[-5] in recount:
			recount[file[-5]] += 1
		if file [-6] in recount:
			recount[file[-6]] += 1
		if file[-5] not in recount:
			recount.setdefault(file[-5], 1)
		if file[-6] not in recount:
			recount.setdefault(file[-6], 1)
	recount_sorted = {k: v for k, v in sorted(recount.items(), key=lambda item: item[1], reverse=True)} # Sorts the dictionary by its values, which represent the most frequent chain in pdb files from input_directory.
	chains = list(recount_sorted.keys())
	core_chain = chains[0] # Select the most frequent chain.
	for file in list_of_pdb_files:
		search_pattern = re.compile(core_chain) # Saves the chain id as a search pattern.
		m = search_pattern.search(file, 5, 7) # Search the most frequent chain through all pdb files
		if m:
			core_chain_file = file # Select the first pdb file contaning the search pattern as core chain file.
			break

	return core_chain_file


def generate_output_file(final_model,out_name):
	""" This function takes as input both the final model created with the building algorithm and the output filename given by the user (if not defined, is macrocomplex by default). Eventually, it returns the file saved in either ".pdb" or ".mmcif" format. """
	out_name = str(out_name.strip())
	# If the output file is too big, we save it in ".mmcif" format
	if len(list(final_model[0].get_atoms())) > 99999 or len(list(final_model[0].get_chains())) > 62:
		mmcif_IO = MMCIFIO()
		mmcif_IO.set_structure(final_model[0])
		mmcif_IO.save(out_name + ".cif")
	# Otherwise, save it ".pdb" format
	else:
		pdb_IO = PDBIO()
		pdb_IO.set_structure(final_model[0])
		pdb_IO.save(out_name + ".pdb")

def superimpose_structures (core_model, test_model, RMSD_threshold):
	""" This function performs the superimpostion of the atoms of every combination between the core and the test chains and returns a dictionary with the superimposed chains as keys and the superimposition objects as values, the best RMSD of those combinations, and a boolean (will be True if there has been at least one superimposition).

	Key arguments:
	core_model -- first model of the core structure object.
	test_model -- first model of the test structure object.
	RMSD_threshold -- root-mean-square deviation threshold for the superimposition.
	"""
	# Declares variables
	best_RMSD = 0
	previous_RMSD = True
	superimposed_chains = False
	superimpositions = {}

	for core_chain in core_model.get_chains(): # Iterates through all the chains of the core model
		#### Obtains the molecule type and the atoms of core chain ####
	    core_molecule_type = get_molecule_type(core_chain)
	    if core_molecule_type == "Protein":
	        core_ca_atoms = get_backbone_atoms_protein(core_chain)
	        core_atoms = (core_ca_atoms, core_molecule_type)
	    else:
	        core_c4_atoms = get_backbone_atoms_nucleicacids(core_chain)
	        core_atoms = (core_c4_atoms, core_molecule_type)

	    for test_chain in test_model.get_chains(): # Iterates through all the chains of the test model
			#### Obtains the molecule type and the atoms of test chain ####
		    test_molecule_type = get_molecule_type(test_chain)
		    if test_molecule_type == "Protein":
			    test_ca_atoms = get_backbone_atoms_protein(test_chain)
			    test_atoms = (test_ca_atoms, test_molecule_type)
		    else:
			    test_c4_atoms = get_backbone_atoms_nucleicacids(test_chain)
			    test_atoms =(test_c4_atoms, test_molecule_type)

		    if core_atoms[1] != test_atoms[1]: # Checks out if both chains are the same type of molecule
			    pass # poner verbose

		    elif len(core_atoms[0]) != len(test_atoms[0]):  # Checks out if both chains have the same number of atoms
			    pass # poner verbose
		    else: # If both chains are the same chain, then is possible to superimpose them.
			    superimposition = Superimposer() # Creates the Superimposer object
			    superimposition.set_atoms(core_atoms[0], test_atoms[0]) # Performs the superimposition of core and test atoms.
			    RMSD_value = superimposition.rms # Saves the RMSD value of the superimposition
			    if RMSD_value > RMSD_threshold: # If the RMSD of the superimposition is bigger than the rmsd threshold, skip this superimposition.
				    continue
			    if previous_RMSD is True or RMSD_value < previous_RMSD: # If RMSD is lower than the threshold, this condition will be checked. If it is true, stores the current RMSD value.
				    previous_RMSD = RMSD_value
				    best_RMSD = RMSD_value

			    superimpositions[(core_chain.get_id(),test_chain.get_id())] = superimposition # Stores in superimpositions dictionary the chain IDs as keys, and its respective result of the superimposition as value.
			    superimposed_chains = True
	if superimposed_chains is True: # If at least has been one superimposition enterns the if condition.
	    superimpositions = sorted(superimpositions.items(), key = lambda x: x[1].rms) # Sorts the dictionary according to RMSD value of the superimpositions.
	return (superimpositions, best_RMSD, superimposed_chains)



def BioBuilder(core_structure, files_list, num_iteration, stop_counter, id_counter, extra_arguments):
	""" This function plays on the superimposition function to construct recursively the macrocomplex biomolecule, and eventually returns it as a structure object when a given condition is met.

	Keyword arguments:
	core_structure -- core molecule as a structure object.
	files_list -- list of input ".pdb" files.
	num_iteration -- number that keeps track of the number of iterations.
	stop_counter -- condition to stop the recursive function.
	id_counter -- support for ID nomenclature discrepances.
	extra_arguments -- Argparse object containing all the arguments entered by the user in the cmd.
			stoichiometry -- number of chains the complex must have.
			RMSD_threshold -- root-mean-square deviation threshold for the superimposition.
			number_clashes -- maximum number of clashes allowed to consider that one chain has not been added yet to the macrocomplex.
			input_directory -- input directory where all pdb files are stored.
			verbose -- to show the progress of the function.
	"""
	# Saves the extra arguments into different internal variable names.
	args = extra_arguments
	pdb_files = files_list
	sto = extra_arguments.stoichiometry
	RMSD_threshold = extra_arguments.RMSD_threshold
	number_clashes = extra_arguments.number_clashes
	input_directory = extra_arguments.input_directory
	verbose = extra_arguments.verbose

	current_num_core_chain = len(core_structure[0])
	num_iteration += 1

	##### Conditions to stop the recursive function ######
	if current_num_core_chain == sto or stop_counter > len(pdb_files):
		return core_structure # the program finishes.
	# if sto != None:
	# 	if current_num_core_chain == sto: # Stops the iteration process if the current number of chains is equal to stoichiometry.
	# 		return core_structure
	# else:
	# 	if stop_counter > len(pdb_files): # or the iterations arrives to the number of pdbs available.
	# 		return core_structure


	test_file = pdb_files[0] # Selects the test file.
	test_path = input_directory + "/" + test_file
	test_structure = obtain_structure("test_file", test_path)
	test_model = test_structure[0]


	superimpositions, best_RMSD, superimposed_chains = superimpose_structures(core_structure[0], test_structure[0], RMSD_threshold) # Assigns the results of superimpose function into different internal variables.


	if superimposed_chains is False or best_RMSD > RMSD_threshold:  # If the superimposition does not reach the conditions to be valid, or there is no superimposition at all, the test file is send to the bottom of the pdb list.
		removed_file = pdb_files.pop(0)
		pdb_files.append(removed_file)
		num_iteration += 1
		stop_counter += 1
		return BioBuilder(core_structure = core_structure, files_list = pdb_files, num_iteration = num_iteration, stop_counter = stop_counter, id_counter = id_counter, extra_arguments = args) # Calls the recursive function.

	else: # If there is a valid superimposition.
		for chains, superimp_obj in superimpositions: # Iterates through all the chains and superimposition objects in superimposition dictionary.
			superimp_obj.apply(test_model.get_atoms()) # Applies the rotation and translation matrix to the atoms of the test model
			adding_chain = []
			for chain in test_model.get_chains(): # Iterates through all the chains in test model.
				if chain.get_id() != chains[0]: # Selects only the chain from the test model which has not been superimposed with the core chain (potential chain to add to the core model).
					test_molecule = get_molecule_type(chain) # Obtains the molecule type in order to get the backbone atoms.
					adding_chain = chain # Stores the potential chain to add.
			chain_already_added = False
			# Gets the backbone atoms #
			if test_molecule == "Protein":
				test_atoms = get_backbone_atoms_protein(adding_chain)
			else:
			 	test_atoms = get_backbone_atoms_nucleicacids(adding_chain)

			for chain in core_structure[0].get_chains(): # Iterates through the chains in core strucuture.
				core_molecule = get_molecule_type(chain)
				if core_molecule == "Protein":
					core_atoms = get_backbone_atoms_protein(chain)
				elif core_molecule == "DNA" or core_molecule == "RNA":
				 	core_atoms = get_backbone_atoms_nucleicacids(chain)

				clash_search = NeighborSearch(core_atoms) # Creates a NeighborSearch object with the atoms of the core chain.
				atom_clash = []
				for atom in test_atoms: # Iterates through all the atoms in the test chain (potential chain to add to the core model).
					list_atoms_clash = clash_search.search(atom.coord, 5) # For each one of the test atoms, looks for clashes with any of the core atoms, and returns a list with those cases that are not under the threshold.
					if len(list_atoms_clash) != 0: # If there are clashes in the list_atoms_clash, adds this list to the list of clashes created.
						atom_clash.extend(list_atoms_clash)
				if len(atom_clash) > number_clashes: # Checks if the number of clashes in the list of clashes is above the threshold.
					chain_already_added = True # Then, the potential chain to add is consedered as added, so we skip that chain.
					if args.verbose:
						sys.stderr.write("Chain %s of the test file %s clashes with chain %s in the complex.\n" %(adding_chain.id, test_file, chain.id))
					break
				elif len(atom_clash) <= number_clashes:
					if args.verbose:
						sys.stderr.write("Chain %s of the test file %s is not clashing with any chain in the complex.\n" %(adding_chain.id, test_file))
					continue
			if chain_already_added is False: # If the potential chain to add is not already present.
				id_counter += 1
				for chain in core_structure[0].get_chains(): # Checks if the ID of the adding_chain is already taken by any of the core chains, in that case, change the adding_chain ID.
					if adding_chain.id == chain.id:
						old_id_chain = adding_chain.id
						# adding_chain.id = id_counter
						adding_chain.id = old_id_chain + str(id_counter)
						if args.verbose:
							sys.stderr.write("ID of the current chain has changed from %s to %s.\n" %(old_id_chain, adding_chain.id))
				core_structure[0].add(adding_chain) # Adds the potential chain to add to the core structure.
				if args.verbose:
					sys.stderr.write("The chain %s has now been added to the macrocomplex.\n" %(adding_chain.id))
				removed_file = pdb_files.pop(0) # Test file is sent to the end of the list.
				pdb_files.append(removed_file)
				num_iteration += 1
				stop_counter = 0 # Restarts the stop_counter variable since there has been a chain addition.
				return BioBuilder(core_structure = core_structure, files_list = pdb_files, num_iteration = num_iteration, stop_counter = stop_counter, id_counter = id_counter, extra_arguments = args) #Calls the recursive function again.
	## No valid superimposition ##
	removed_file = pdb_files.pop(0) # Test file is sent to the end of the list, in order to get a new one.
	pdb_files.append(removed_file)
	num_iteration += 1
	stop_counter += 1
	return BioBuilder(core_structure = core_structure, files_list = pdb_files, num_iteration = num_iteration, stop_counter = stop_counter, id_counter = id_counter, extra_arguments = args) # Calls the recursive function again.
