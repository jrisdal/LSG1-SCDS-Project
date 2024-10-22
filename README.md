This set of programs was created by Jordan Risdal for completion of the Scientific Computation and Data Sciences Certificate at the University of Texas at Austin utilizing data generated in Dr. Arlen Johnson's lab.

The input data in this case is sam files produced from bowtie2. Sam files are from next gen sequencing data produced through PCR mutagenesis with reads of 75nt. Programs should be run in the following order in order to count mutational variation, clean output, calculate frequencies of amino acid substitutions, and obtain fitness scores for each substitution and tolerance scores for each residue

mut_count.py > output_cleaner.py > AA_freq.py > fitness.py & tolerance.py

heatmap.R is for visualization of fitness scores in a heatmap, and LSG1 pymol.py is for mapping tolerance scores onto preexisting structure in pymol
