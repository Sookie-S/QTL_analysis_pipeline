For more information, please contact soukaina.timouma@manchester.ac.uk


############################################
0- There are initially 2 folders, Files and Scripts, in the root directory (for example here, the root folder is 'Rachel/Sookie/') and a Results folder that will be created by the pipeline:
############################################

- a "Files" folder:

	* It contains the folders for the conditions tested (for example "ethanol", "fluconazole", etc), that have each two sub-directories for each parental species that contains the LOD measured in each chromosome (.out files).
	
        * There is also a 'Genomes' folder that contains a folder called GFF_files where you need to store the annotation of parent 1 and parent 2 (GFF files).

- a Scripts folder:

        * contains all the scripts for this analysis

- a "Results" folder:

	* the folder 'QTL_regions' contains sub-folders for each condition tested (for example "ethanol", "fluconazole", etc), and it contains the QTL regions detected (.csv files)

	* the folder 'Gene_overlap' contains sub-folders for each condition tested (for example "ethanol", "fluconazole", etc), and it gives as information the genes that overlap with each QTL regions, and the information if they are within the LOD support interval.


############################################
1- Make sure the input files (.out files) for parent1 and parent2 are in two separated folders, that contains "parent1" or "parent2" in their names. Make sure that there is only the two folders containing the .out files, the folder for parent1 and the for parent2.
############################################

For example, for the hybrid S. cerevisiae/S. kudriavzevii, is parent1 name is "scer" and parent2 name is "skud", the names of the folders can be:

- Updated_Results_lowtemp_skud_250222
- Updated_Results_lowtemp_scer_250222
or
- MPresults_scer
- MPresults_skud
or
- WhateverBLAscer
- WhateverBLAskud


############################################
2- Make sure that the folders containing the .out files are named appropriately: "SUFFIX_CHR.out".
############################################

The suffix "SUFFIX" can be named anything as long as it doesn't contain an "_"
The prefix "CHR" corresponds to the chromosome name as written in the GFF file.

for example, for the hybrid S. cerevisiae/S. kudriavzevii:
The Scer chromosomes "CHR" should be named chrI, chrII, chrIII etc, as it is the way they are named the GFF files.
The Skud chromosomes "CHR" should be named Skud_1, Skud_2, Skud_3 etc, as it is the way they are named the GFF files.

For example:
  Valid naming of the output files:
   - scerLT_chrI.out
   - skidLT_Skud_1.out
   - whatever_chrI.out
   - whatever_Skud_1.out
  Invalid naming (TO NOT DO, it's impertative there is no "_" in the SUFFIX):
   - scer_LT_chrI.out
   - skidLT_lowtemp_bla_bli_blu_Skud_1.out
   - scer_chr1.out
   - skidLT_Skud1.out


############################################
3- launching the pipeline (from the script folder directory)
############################################


To launch the pipeline, please specify all the following arguments:

./pipeline.sh [LOD_threshold] [LOD_support_interval] [conditionName] [Name_Parent1] [Name_parent2] [ChromosomePrefix_Parent1] [GFF_file_Name_Parent1] [ChromosomePrefix_Parent2] [GFF_file_Name_Parent2]


for example:
- LOD_threshold = 5
- LOD_support_interval = 1
- conditionName = test
- Name_Parent1 = scer
- Name_parent2 = skud
- ChromosomePrefix_Parent1 = chr
- GFF_file_Name_Parent1 = YPS128.all_feature.gff
- ChromosomePrefix_Parent2 = Skud
- GFF_file_Name_Parent2 = skud_updated.gff

That would be:
./pipeline.sh 5 1 test scer skud chr YPS128.all_feature.gff Skud skud_updated.gff

