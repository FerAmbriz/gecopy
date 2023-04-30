#!/bin/bash
#Ask the user for the directory where the tumor bams are stored
echo "Tumor bams directory:" $1
#Same sequence is repeated for every variable necessary to run CNVkit
#Normal files
echo "Normal bams directory:"$2
#Targets bed file
echo "Path to your targets bed file:"$3
#RefFlat file
echo "Path to your refFlat.txt file (provided by the CNVkit package):" $4
#Fasta reference file
echo "Path to your fasta reference genome file:" $5
#Mappable access file
echo "Path to your mappable access file:" $6
#CNVkit reference name and path
echo "Where do you want me to store the CNVkit reference?:" $7
echo "How should I name the reference?" $8
#Output directory file
echo "Where do you want me to store the CNVkit results files?:" $9
# ******** Parameters specific to EXOMEDEPTH ************
# Números de más de un dígito van entre {} 
echo "Path to the adapted bedfile to be used for the EXOME DEPTH analysis: " ${10}
echo "Path to the directory where you want to store EXOMEDEPTH results..." ${11}
echo "Path to the adapted bedfile to be used for the Panelcn.mops analysis:" ${12}
echo "List the genes of interest for this analysis (separated by a ','):" ${13}
echo "Name for Panelcn.mops results file:" ${14}
echo "Path to the directory where you want to store Panelcn.mops results..." ${15}
echo "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" "${12}" "${13}" "${14}" "${15}"
# --------------------------------------------------------------------------
# Then we run CNVkit with the variables inserted by the user
echo "----------------------------------------------------------------------"
echo "Initializing CNVkit analysis..."
cnvkit.py batch "$1"/*.bam --normal "$2"/*.bam \
    --targets "$3" --annotate "$4" \
    --fasta "$5" --access "$6" \
    --output-reference "$7"/"$8".cnn  --output-dir "$9" \
    --diagram --scatter
echo "CNVkit general pipeline... DONE"
# ---------------------------------------------------------------------------
# We proceed to the calling process for each sample results
echo "-----------------------------------------------------------------------"
echo "Initializing CNVkit calling pipeline"
for entry in "$9"/*.cns
do
  cnvkit.py call "$entry" -o "$entry".call
done
echo "DONE"
# ---------------------------------------------------------------------------
# Finally we perform the genemetrics pipeline for CNVkit
echo "-----------------------------------------------------------------------"
echo "Initializing CNVkit genemetrics analysis"
for entry2 in "$9"/*.cnr
do
	cnvkit.py genemetrics "$entry2" > "$entry2".genemetrics.txt
done
echo "DONE"
# ---------------------------------------------------------------------------
echo "***********************************************************************"
echo "CNVkit analysis... DONE"
echo "***********************************************************************"
# ---------------------------------------------------------------------------
echo "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" "${12}" "${13}" "${14}" "${15}"
echo "Initializing CNV analysis with EXOME DEPTH..."
Rscript Exome_Depth_bash.R "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}"
echo "Analysis of CNV with EXOME DEPTH... DONE"
# ---------------------------------------------------------------------------
echo "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" "${12}" "${13}" "${14}" "${15}"
echo "Initializing CNV analysis with Panelcn.mops"
Rscript panelcn_mops_bash.R "$1" "$2" "$3" "$4" "$5" "$6" "$7" "$8" "$9" "${10}" "${11}" "${12}" "${13}" "${14}" "${15}"

echo "Analysis of CNV with Panelcn.mops ... DONE"
echo "CNV analysis... COMPLETED"

