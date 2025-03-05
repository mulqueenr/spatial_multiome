# spatial_multiome
 Nextflow processing of curio+multiome data

## Set up on Geo
```bash
cd /volumes/USR2/Ryan/projects/spatial_wgs
mkdir -p ./tools

git clone https://github.com/mulqueenr/spatial_multiome.git

example run
source activate #to use more recent version of java

#first need to make the output dir and the log directory for bcl-convert
flowcellDir = "/volumes/seq/flowcells/MDA/nextseq2000/2025/250227_RM_CurioWGS_scalemet"
outdir = "/volumes/USR2/Ryan/projects/spatial_wgs/250129_First_Experiment"
mkdir -p ${outdir}
mkdir -p ${outdir}/logs

cd /volumes/USR2/Ryan/projects/10x_MET #move to project directory
git clone https://github.com/mulqueenr/scmet_nf_processing #pull github repo

echo """[Settings],
CreateFastqForIndexReads,1
OverrideCycles,Y50;I8N2;U24;Y47
[Data],
Sample_ID,index
dcis41t,ACGAGTAG
dcis41t,CAATCCCT
dcis41t,GTCCAGGC
dcis41t,TGTGTATA""" > ${outdir}/DNA_SampleSheet.csv

#sequencing_cycles in quotes to avoid newline char

nextflow ./tools/spatial_multiome/spatial_processing.groovy \
-with-report \
--flowcellDir ${flowcellDir} \
--sequencing_cycles="Y50;I8N2;N8I16;Y47" \ 
--outname 250129_spatialdna \
--outdir ${outdir} \
--resume
```