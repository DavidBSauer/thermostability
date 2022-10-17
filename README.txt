# Thermostability
Script of two bioinformatic methods to identify a limited set of amino acids or positions that likely underlie thermostability from growth temperature data and protein multiple sequence alignments. 

See: Sauer, Karpowich, Song & Wang. Rapid Bioinformatic Identification of Thermostabilizing Mutations. Biophysical Journal (2015) https://doi.org/10.1016/j.bpj.2015.07.026

## Requirements
The script requires Biopython, NumPy, MatPlotLib, and matplotlib-venn be installed. 

http://biopython.org/wiki/Download
http://www.scipy.org/scipylib/download.html
http://matplotlib.org/downloads.html
https://pypi.python.org/pypi/matplotlib-venn

## Running the script
Start the script by running in the terminal. After starting, the script will ask the analysis type, if appropriate, and appropriate parameters. The script will then run automatically, with ETAs printed for long steps.

The input type can be either species assignment only (SAO), aligned (ALN), or species and temps assigned (STA).
	SAO will take the multiple sequence alignment (MSA) file, assign species, and only retain those sequences with species information (helpful to get a sense of the OGT range quickly, or to more efficiently realign only the useful sequences). After starting everything will run automatically.
	run as: python automated_whole.py W.fa X.txt Y SAO
	where: 	W MSA in fasta format
		X is a CSV file of species (first column) and OGTs (second column)
		Y is the run name, used for creating sub-directories 

	ALN will take an input alignment file and look up species for each sequene, assign temperatures, calculate identities, and run the analysis. 
	run as: python automated_whole.py W.fa X.txt Y ALN
	where: 	W input in fasta format
		X is a csv file of species (first column) and OGTs (second column)
		Y is the run name, used for creating sub-directories 

	STA will take an input alignment file and identity matrix and run the analysis. This is particularly helpful in repeating an analysis without reassigning species and recalculating identities. 
	run as: python automated_whole.py W.fa X.txt Y STA
	where: 	W input in fasta format
		X is a identity matrix (Clustal format)
		Y is the run name, used for creating sub-directories 

## Result files and figures
There will be one folder created using the given run name. Inside of that will be basic data organization files and, if appropriate, two subdirectories for the “global” and “pairwise” methods. 

In the “Y” folder will be (where Y is the run name):
	Y_ident_matrix.txt - an identity matrix for the sequences, in clustal format.
	Y_raw_seq.fa - the alignment file with assigned sequences, in fasta format.
	Y_Xed.fa - the alignment file, with ambigious amino acids removed, in fasta format.
	Y.log - a log file with notes on each step of the analysis.
	temps_assigned.png - a histogram of sequences with assigned temperatures. 
	Y_species_temp_data_used_venn.png - a Venn diagram of species-OGTs used out of all OGTs provided.
	Y_sequences_assigned_of_input_venn.png - a Venn diagram of sequences assigned an OGT out of all the sequences provided.
	Y_

In the “Y/pairwise” folder will be (where Y is the run name):
	temp_dependent_differences.csv - a CSV file of the positional frequency of amino acid differences, by position, for those sequence pairs that are sufficiently identical and with sufficiently different OGTs.
	references_differences.csv - a CSV file of the positional frequency of amino acid differences, by position, for those sequence pairs that are sufficiently identical.
	temp_dependent_minus_ref_differences.csv - a CSV file of the temperature depended differences only. (temp_dependent_differences - reference_differences, by position).
	gap_frequency.csv - a CSV file of the positional gap frequency of the MSA.
	temp_dependent_minus_ref_differences.csv - a CSV file of the temperature depended differences only. Positions with a gap frequency greater than the provided threshold are set to zero.
	higher_residues.csv - a CSV file of the amino acids of the higher temperature amino acid differences, by position. 
	lower_residues.csv - a CSV file of the amino acids of the lower temperature amino acid differences, by position. 

In the “Y/pairwise/figures” folder will be (where Y is the run name):
	temp_dependent_differences.png - Bar chart of the positional frequency of amino acid differences, by position, for those sequence pairs that are sufficiently identical and with sufficiently different OGTs.
	temp_dependent_differences.png - Bar chart of the positional frequency of amino acid differences, by position, for those sequence pairs that are sufficiently identical and with sufficiently different OGTs.
	references_differences.png - Bar chart file of the positional frequency of amino acid differences, by position, for those sequence pairs that are sufficiently identical.
	temp_dependent_minus_ref_differences.png - Bar chart of the temperature depended differences only. (temp_dependent_differences - reference_differences, by position).

In the “Y/pairwise/figures/histograms” folder will be (where Y is the run name) a series of files named “pos_x_histogram.png” (where X is the column of the MSA) of the higher (red) and lower (blue) temperature amino acids frequencies. 

In the “Y/global” folder will be (where Y is the run name):
	Y_non_redundant.fa - a MSA of the sequences used for the global analysis (after removing redundant sequences, based on the provided identity threshold).
	hyperthermophiles.fa, thermophiles.fa, mesophiles.fa, psychrophiles.fa - MSAs of sequences in each OGT range (if any).
	hyper.csv, thermo.csv, meso.csv, psych.csv - A CSV file of amino acid counts, by position, for the OGT ranges of hyperthermophiles, thermophiles, mesophiles, and psychrophiles respectively (if any).
	hyper_norm.csv, thermo_norm.csv, meso_norm.csv, psych_norm.csv - A CSV file of amino acid frequencies, by position, for the OGT ranges of hyperthermophiles, thermophiles, mesophiles, and psychrophiles respectively (if any).
	hyper+thermo_norm.csv - A CSV file of amino acid frequencies, by position, for thermophiles and hyperthermophiles combined.
	hyper-meso.csv - A CSV file of difference in the amino acid frequencies, by position, between hyperthermophiles and mesophiles.
	hyper+thermo-meso.csv - A CSV file of difference in the amino acid frequencies, by position, for thermophiles and hyperthermophiles combined versus mesophiles.

In the “Y/global/figures” folder will be (where Y is the run name):
	global_temps_used_histogram.png - A histogram of the temperatures in the non-redundant MSA (Y_non_redundant.fa) used in this analysis.
	redundancy_venn.png - Venn diagram of the sequences used relative to the total sequences provided from the OGT assigned MSA.
	hyper_norm.png, thermo_norm.png, meso_norm.png, psych_norm.png - Heat map of amino acid frequency, by position, for hyperthermophiles, thermophiles, mesophiles, and psychrophiles respectively (if any).
	hyper_thermo_norm.png - Heat map of amino acid frequency, by position, for hyperthermophiles and thermophiles combined.
	hyper_meso_diff.png - Difference heat map of amino acid frequencies, by position, between hyperthermophiles and mesophiles.
	hyper_thermo_meso_diff.png - Difference heat map of the amino acid frequencies, by position, for thermophiles and hyperthermophiles combined versus mesophiles.

In the “Y/global/figures/histograms” folder will be (where Y is the run name)a series of files named “pos_x_histogram.png” (where X is the column of the MSA) of the hyperthermophile and thermophile combined (red) and mesophile amino acid frequencies.

## General Notes
The analysis is absolutely dependent on alignment quality. Make sure the seuqences are well aligned before the analysis. Pathologies are easy to note if systematic difference are seen in the global heatmaps. 
The calculation time scales with alignment length, removing poorly aligned portions will significantly shorten the calculation time. 
The pairwise calculation takes significantly longer than the global analysis, but is calculated first when run together. If earlier results are desired, run the global analysis seperately. 
	
