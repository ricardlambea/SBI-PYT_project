# Example Package

This is a simple example package. You can use
[Github-flavored Markdown](https://guides.github.com/features/mastering-markdown/)
to write your content.

First of all it is important to state that we have created a directory "examples" in which we have all the folders that are used as input directories for our program, and that the command is always executed from within that "examples" directory.


EXAMPLE 1:
The first example is the protein 1gzx, which is the oxy T state haemoglobin from Homo sapiens. This protein hetero 4-mer (A2B2).
To run our program we execute the following command:

  python3 biobuilder_core.py -i /1gzx -o out_dir -of 1gzx -sto 4 -v

Note that we are executing our program with python3 at the beginning, and that is because the program has not been installed, but if it has been installed the program, "python3" would not be necessary. That works for all the examples.

About the arguments used after the script name, there are the mandatory input (-i) and output (-o) directories, we also called the optional output filename (-of) argument, the stoichiometry (-sto) argument with a value of 4 as it is the number of chains present in the complex, and finally the verbose argument (-v) to redirect to the standard error channel the verbose statements.

Using "time" at the beginning of the command ("time python3 biobuilder...") when running the shell we can see the total amount of time the process lasted. In this example it took 0.810 seconds.

Matchmaker 1gzx_original.pdb, chain B (#1) with 1gzx.pdb, chain B (#0), sequence alignment score = 788.5
with these parameters:
	chain pairing: bb
	Needleman-Wunsch using BLOSUM-62
	ss fraction: 0.3
	gap open (HH/SS/other) 18/18/6, extend 1
	ss matrix:  (O, S): -6 (H, O): -6 (H, H): 6 (S, S): 6 (H, S): -9 (O, O): 4
	iteration cutoff: 2
RMSD between 146 pruned atom pairs is 0.000 angstroms; (across all 146 pairs: 0.000)


EXAMPLE 2:
The second example is the protein 3kuy, which corresponds to DNA stretching in the nucleosome, which in turn facilitates alkylation by an intercalating antitumor agent, and it comes from Escherichia coli. This protein is a hetero 8-mer (A2B2C2D2).
To build our model we execute the following command:

  python3 biobuilder_core.py -i /3kuy -o out_dir -of 3kuy -sto 10 -v

In that case we use 10 as the stoichiometry value as there are 8 protein chains and 2 acid nucleic molecules. The computational time this example took was 2.989 seconds. As can be seen, the superimposition is not perfect, but it is quite good.

En est caso nos da 0 A el RMSD pero no me cuadra, la superposicion no es perfecta en chimera.



EXAMPLE 3:
The third example is the protein 5ara, which is a bovine mitochondrial ATP synthase from E. coli BL21(DE3). This protein is a hetero 22-mer (A8B3C3DEFGHIJK).
In order to buil the complex we run on the shell the command:

  python3 biobuilder_core.py -i /5ara -o out_dir -of 5ara -sto 22 -v

The computational time in this case was 13.222 seconds.

Matchmaker 5ara_original.pdb, chain A (#1) with 5ara.pdb, chain A (#0), sequence alignment score = 2557.5
with these parameters:
	chain pairing: bb
	Needleman-Wunsch using BLOSUM-62
	ss fraction: 0.3
	gap open (HH/SS/other) 18/18/6, extend 1
	ss matrix:  (O, S): -6 (H, O): -6 (H, H): 6 (S, S): 6 (H, S): -9 (O, O): 4
	iteration cutoff: 2
RMSD between 509 pruned atom pairs is 0.000 angstroms; (across all 509 pairs: 0.000)



EXAMPLE 4:
The fourth example is the protein 5dn6, an ATP synthase from Paracoccus dentrifricans (strain Pd 1222). It is a hetero 27-mer (A12B3C3DEFGHIJKL).
To build the model we execute the following command in the shell:

  python3 biobuilder_core.py -i /5dn6 -o out_dir -of 5dn6 -sto 27 -v

The computational time was 14.164 seconds.



EXAMPLE 5:
The fifth example is the protein 5oom, a structure of a native assembly intermediate of the human mitochondrial ribosome with unfolded interfacial rRNA. It is a hetero 53-mer  (ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzA).
To build the model we execute the following command:

  python3 biobuilder_core.py -i /5oom -o out_dir -of 5oom -sto 53 -v

The computational time has been 1 minute and 7.917 seconds.



EXAMPLE 6:
The sixth example is the protein 6ezm, which is a imidazoleglycerol-phosphate dehydratase from Saccharomyces cerevisiae. It is a homo 24-mer (A24).
We run the next command on the shell to build it:

  python3 biobuilder_core.py -i /6ezm -o out_dir -of 6ezm -sto 24 -v

The computational time has been 6.389 seconds.



EXAMPLE 7:
The seventh example is the protein 5vox, which is a V-ATPase from S. cerevisiae (strain ATCC 204508/S288c). It is a hetero 33-mer (A8B3C3D3E3F3GHIJKLMNOP).
To build the model we execute the following command:

  python3 biobuilder_core.py -i /5vox -o out_dir -of 5vox -sto 33 -v

The computational time was 22.505 seconds.

Matchmaker 5vox_original.pdb, chain b (#1) with 5vox.pdb, chain b (#0), sequence alignment score = 3058.7
with these parameters:
	chain pairing: bb
	Needleman-Wunsch using BLOSUM-62
	ss fraction: 0.3
	gap open (HH/SS/other) 18/18/6, extend 1
	ss matrix:  (O, S): -6 (H, O): -6 (H, H): 6 (S, S): 6 (H, S): -9 (O, O): 4
	iteration cutoff: 2
RMSD between 634 pruned atom pairs is 0.000 angstroms; (across all 634 pairs: 0.000)
