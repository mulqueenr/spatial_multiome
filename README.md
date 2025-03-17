# spatial_multiome
 Nextflow processing of curio+multiome data
 Running on department server

## Set up on Geo
```bash
cd /home/rmulqueen/projects/spatial_wgs/
mkdir -p ./tools

git clone https://github.com/mulqueenr/spatial_multiome.git ./tools/spatial_multiome

#example run

#first need to make the output dir and the log directory for bcl-convert
DNA_flowcellDir="/home/rmulqueen/projects/spatial_wgs/seq/250227_RM_CurioWGS_scalemet"
RNA_flowcellDir="/home/rmulqueen/projects/spatial_wgs/seq/250220_RM_CuioWGS_RNA"

outdir="/home/rmulqueen/projects/spatial_wgs/data/250129_First_Experiment"
outname="dcis41t"
mkdir -p ${outdir}
mkdir -p ${outdir}/logs

cd /home/rmulqueen/projects/spatial_wgs/ #move to project directory
git clone https://github.com/mulqueenr/spatial_multiome.git ./tools/spatial_multiome #pull git repo

echo """[Settings],
CreateFastqForIndexReads,1
OverrideCycles,Y50;I8N2;N8U16;Y47
[Data],
Sample_ID,index
${outname}_dna,ACGAGTAG
${outname}_dna,CAATCCCT
${outname}_dna,GTCCAGGC
${outname}_dna,TGTGTATA""" > ${outdir}/DNA_SampleSheet.csv

echo """Lane,Sample,Index
*,${outname}_spatial,SI-TT-D6
*,${outname}_rna,SI-TT-D11""" > ${outdir}/RNA_SimpleSampleSheet.csv

#sequencing_cycles in quotes to avoid newline char
#install conda from conda-forge

source activate #just to make sure i'm using a more up-to-date java

nextflow ./tools/spatial_multiome/nextflow_running/spatial_processing.groovy \
-with-report \
--dna_flowcellDir ${DNA_flowcellDir} \
--dna_samplesheet ${outdir}/DNA_SampleSheet.csv \
--rna_flowcellDir ${RNA_flowcellDir} \
--rna_samplesheet ${outdir}/RNA_SimpleSampleSheet.csv \
--outname ${outname} \
--outdir ${outdir}
```