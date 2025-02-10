#Repository of code for "Diversity of B Cell Populations and Ig Repertoire in Human Lungs. J Immunol. 2023 Aug 1;211(3):486-496. doi: 10.4049/jimmunol.2200340. PMID: 37314411; PMCID: PMC10352589."

Cha_Pipeline_Code.R
	#Pipeline code that calculates diversity ratios (DR) using commands from the developing Cha package.
	##Commands use Fumi version of code that has been pushed onto GitHub
	##Under Fumi version, use of iNEXT estimates in DR can be toggled on or off.
	##HillRatio function code to loop across samples was lost during data transfer
	##Function codes in Sequencing_Function_Codes
	##Read input files from Cloanalyst

Sequencing_Function_Codes.R
	#All function codes used for analyzing PartitionAnalysis files

Diversity_Ratio_V1
	#Test output file from Cha_Pipeline_Code.R
	##Contains DRs from each patient sample using iNEXT estimates.

Sequencing_Analysis_Code.R
	#Code for manually checking for quality control in samples
	##Checks for clone degeneracy, top QID, Entrophy, Jaccard index, Simpson index, INEXT, delta of clone divergence.
