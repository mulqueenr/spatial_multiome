//Nextflow pipeline for processing Navin lab spatial multiome//
// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
// DNA
params.dna_flowcellDir = "/home/rmulqueen/projects/spatial_wgs/seq/250227_RM_CurioWGS_scalemet" //Sequencing run flowcell dir
params.dna_samplesheet = "DNA_SimpleSampleSheet.csv"
params.dna_bases_mask = "Y50,I8N*,Y24,Y47"
params.dna_flowcell = "AAGHWKTM5"

// RNA
params.rna_flowcellDir = "/home/rmulqueen/projects/spatial_wgs/seq/250220_RM_CuioWGS_RNA" //Sequencing run flowcell dir
params.rna_samplesheet = "RNA_SimpleSampleSheet.csv"
params.rna_bases_mask = "Y28,I10,I10,Y90"
params.rna_flowcell = "AAGHWVCM5"
params.spatial_barcode = "/home/rmulqueen/tools/curiotrekker-v1.1.0/U0028_003_BeadBarcodes.txt"

//REF
params.ref="/home/rmulqueen/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
params.src="/home/rmulqueen/projects/spatial_wgs/tools/spatial_multiome/src"
params.cellranger_arc="/home/rmulqueen/tools/cellranger-arc-2.0.2/cellranger-arc"
params.curio_trekker="/home/rmulqueen/tools/curiotrekker-v1.1.0/"

params.max_cpus="200"

//output
params.outname = "250129_spatialdna"
params.outdir = "/home/rmulqueen/projects/spatial_wgs/data/250129_First_Experiment2"
params.date = "20250226"
//library parameters
params.cell_try="5000" //Based on expected cell count from library generation
params.samplesheet="/home/rmulqueen/projects/spatial_wgs/data/250129_First_Experiment/DNA_SampleSheet.csv" //Based on expected cell count from library generation

log.info """

		================================================
		        Spatial Multiome Pipeline v1.1
		================================================
		____________________DNA_________________________
		DNA Flowcell Dir : ${params.dna_flowcellDir}
		DNA Sample Sheet : ${params.dna_samplesheet}

		____________________RNA_________________________
		RNA Flowcell Dir : ${params.rna_flowcellDir}
		RNA Sample Sheet : ${params.rna_samplesheet}

		____________________PARAMS______________________
		Output Directory : ${params.outdir}
		Output Prefix : ${params.outname}
		NF Working Dir : ${workflow.launchDir}
		Cellranger ARC install : ${params.cellranger_arc}

		Max cpus : ${params.max_cpus}
		================================================

""".stripIndent()

// BCL TO FASTQ PIPELINE FOR GENERATING SINGLE-CELL FASTQs
process DNA_CELLRANGER_MKFASTQ { 
	//Generate Undetermined Fastq Files from BCL Files.
	//bcl-convert requires write access to "/var/logs/bcl-convert", so we just bind a dummy one if we add a sif
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/dna_fq", mode: 'copy', overwrite: true

	input:
		path(dna_flowcellDir)
		path(dna_samplesheet)
	output:
		path("dna_fq/${params.dna_flowcell}/${params.outname}_wgs/*fastq.gz"), emit: dna_fq
    script:
		"""
        ${params.cellranger_arc} \\
        mkfastq --id=${params.outname}_dna \\
        --run=${dna_flowcellDir} \\
        --use-bases-mask=${params.dna_bases_mask} \\
        --output-dir=\${PWD}/dna_fq \\
        --samplesheet=${dna_samplesheet}
		"""
}


process RNA_CELLRANGER_MKFASTQ {
	//Run cellranger on RNA samples, this works straight from bcl files.
	// Run mkfastq then run count
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/rna_fq", mode: 'copy', overwrite: true

	input:
		path(rna_flowcellDir)
		path(rna_samplesheet)
	output:
		path("rna_fq/${params.rna_flowcell}/*spatial*fastq.gz"), emit: spatial_fq
		path("rna_fq/${params.rna_flowcell}/*gex*fastq.gz"), emit: gex_fq

    script:
		"""
        ${params.cellranger_arc} \\
        mkfastq --id=${params.outname}_rna \\
        --run=${rna_flowcellDir} \\
        --use-bases-mask=${params.rna_bases_mask} \\
        --output-dir=\${PWD}/rna_fq \\
        --samplesheet=${rna_samplesheet}
		"""

}

process CELLRANGER_COUNT {
	//Run cellranger on DNA samples, to generate GEM-indexed bam file.
	//This throws an error when done, not sure why
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/arc_cellranger_out", mode: 'copy', overwrite: true

	input:
		path(dna_fq), stageAs: "dna_fq/"
		path(gex_fq), stageAs: "gex_fq/"

	output:
		path("${params.outname}/outs/"), emit: multiome_outdir

    script:
		"""
        echo "fastqs,sample,library_type" > libraries.csv
        echo "\${PWD}/gex_fq/,${params.outname}_gex,Gene Expression" >> libraries.csv
        echo "\${PWD}/dna_fq/,${params.outname}_wgs,Chromatin Accessibility" >> libraries.csv

        ${params.cellranger_arc} count \\
		--id=${params.outname} \\
        --reference=${params.ref} \\
        --libraries=libraries.csv \\
        --localcores=${params.max_cpus} \\
        --localmem=1000
		"""

}

process DNA_SPLIT_BAM {
	//Use timoast sinto to split bam by CB tag
	//conda install sinto
	//cpu limit hard set because this causes a lot of i/o
	cpus 50
	publishDir "${params.outdir}/dna/sc_dna_bam", mode: 'copy', overwrite: true 

	input:
		path(multiome_outdir)
	output:
		path("./sc_dna_bam/*bam"), emit: bam

    script:
		"""
		#make splitting barcode list from whitelist
		zcat ${params.outname}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | awk -F, 'OFS="\t" {print \$1,\$1}' > cell_id.tsv

		#split to chunks of 500 cells for i/o purposes
		split -l 500 --numeric-suffixes cell_id.tsv cell_id.split.

		#run cell splitting for each 500 chunk
		for i in cell_id.split* ; do
		sinto filterbarcodes --bam ${params.outname}/outs/atac_possorted_bam.bam --cells \$i -p ${task.cpus} --barcodetag "CB" --outdir ./sc_dna_bam ;
		done
		"""
}


process DNA_PROJECT_COMPLEXITY {
	//Use picard tools to project library complexity
	maxForks 200
	publishDir "${params.outdir}/dna/sc_bam_dedup", mode: 'copy', overwrite: true, pattern: "*bbrd.bam"
	publishDir "${params.outdir}/reports/dna/complexity", mode: 'copy', overwrite: true , pattern: "*.projected_metrics.txt"
	publishDir "${params.outdir}/reports/dna/rmdup", mode: 'copy', overwrite: true , pattern: "*.markdup.log"

	input:
		path(bam)

	output:
		path("*bbrd.bam"), emit: bam_rmdup
		path("*.projected_metrics.txt"), emit: complexity_metrics
		path("*.markdup.log"), emit: rmdup_metrics

    script:
		"""
		cellid=${bam.simpleName}
		#output bam, remove duplicates
		samtools sort -m 10G -n $bam | \\
		samtools fixmate -p -m - - | \\
		samtools sort -m 10G -O BAM -o \${cellid}.tmp.bam

		#output of removal of duplicate reads
		samtools markdup --mode t -r -S -s -f \${cellid}.markdup.log \${cellid}.tmp.bam \${cellid}.bbrd.bam

		#generate library complexity based on 10% downsample rates
		#count unique chr:start sites
		for i in \$(seq 0.1 0.1 1.0); do
		uniq_count=\$(samtools view -F 3332 -s \$i \${cellid}.tmp.bam \\
		| awk 'OFS="\\t"{print \$3,\$4}' \\
		| sort \\
		| uniq -c \\
		| wc -l)
		total_count=\$(samtools view -F 3332 -s \$i \${cellid}.tmp.bam | wc -l)
		echo "\${cellid},\${i},\${total_count},\${uniq_count}"; done > \${cellid}.projected_metrics.txt
		#excluding reads that meet any below conditions:
		#read unmapped (0x4)
		#not primary alignment (0x100)
		#read is PCR or optical duplicate (0x400)
		#supplementary alignment (0x800)
		"""
}

process DNA_COPYKIT {
	// CNV PROFILING 
	//COPYKIT FOR CLONE CALLING BY CNVS
	//Run CopyKit and output list of bam files by clones
	cpus "${params.max_cpus}"
	label 'cnv'
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"
	publishDir "${params.outdir}/dna/cnv_calling", mode: 'copy', pattern: "*{tsv,rds}"
	publishDir "${params.outdir}/plots/cnv", mode: 'copy', pattern: "*pdf"

	input:
		path bam_bbrd_collection
	output:
		path("*.scCNA.rds"), emit: copykit_rds
		path("*.scCNA.tsv"), emit: copykit_tsv
		path("*pdf"), emit: copykit_plots
	script:
	/*
	singularity shell --bind /home/rmulqueen/projects/spatial_wgs/tools/spatial_multiome/src:/src/ ~/singularity/copykit.sif 
	Rscript /src/copykit_cnv_clones.nf.R \
	--input_dir . \
	--output_prefix spatial_dcis41t \
	--task_cpus 300
	*/

		"""
		Rscript /src/copykit_cnv_clones.nf.R \\
		--input_dir . \\
		--output_prefix ${params.outname} \\
		--task_cpus ${task.cpus}
		"""
}

process SPATIAL_CURIO {
	//Run Curio Trekker Pipeline to generate spatial location
	cpus "${params.max_cpus}"
	containerOptions "--bind ${params.src}:/src/,${params.outdir},${params.curio_trekker}:/curio/"
	publishDir "${params.outdir}/spatial", mode: 'copy', overwrite: true

	input:
		path(spatial_barcode)
		path(spatial_fq)
		path(multiome_outdir), stageAs: 'multiome_outdir/*'

	output:
		path("*")

    script:
	"""
	cp \$(realpath ./multiome_outdir/${params.outname}/outs/filtered_feature_bc_matrix) ./filtered_feature_bc_matrix 
	fq1="${params.outname}_spatial_S1_L001_R1_001.fastq.gz"
	fq2="${params.outname}_spatial_S1_L001_R2_001.fastq.gz"

	echo 'sample,sc_sample,experiment_date,barcode_file,fastq_1,fastq_2,sc_outdir,sc_platform,profile,subsample,cores' > samplesheet.trekker.csv
	echo "${params.outname}_rna,${params.outname}_rna,${params.date},${spatial_barcode},\${fq1},\${fq2},\${PWD}/filtered_feature_bc_matrix,TrekkerU_C,singularity,no,${task_cpus}" >> samplesheet.trekker.csv

	bash /curio/nuclei_locater_toplevel.sh \\
	samplesheet.trekker.csv
	"""
}

workflow {
// BCL TO FASTQ PIPELINE FOR SPLITTING FASTQS		
	dna_flowcell_dir = Channel.fromPath(params.dna_flowcellDir)
	dna_samplesheet = Channel.fromPath(params.dna_samplesheet)
	rna_flowcell_dir = Channel.fromPath(params.rna_flowcellDir)
	rna_samplesheet = Channel.fromPath(params.rna_samplesheet)
	spatial_barcode = Channel.fromPath(params.spatial_barcode)

// RUN CELLRANGER ARC PIPELINE
	DNA_CELLRANGER_MKFASTQ(dna_flowcell_dir,dna_samplesheet)
    RNA_CELLRANGER_MKFASTQ(rna_flowcell_dir,rna_samplesheet)
	
	dna_fq = DNA_CELLRANGER_MKFASTQ.out.dna_fq | collect
	
	gex_fq = RNA_CELLRANGER_MKFASTQ.out.gex_fq | collect
	
	spatial_fq = RNA_CELLRANGER_MKFASTQ.out.spatial_fq | collect
    CELLRANGER_COUNT(dna_fq,gex_fq)
	
// DNA PROJECT SINGLE CELL DNA COMPLEXITY
    CELLRANGER_COUNT.out.multiome_outdir \
	| DNA_SPLIT_BAM \
	| DNA_PROJECT_COMPLEXITY

// DNA RUN COPYKIT
	DNA_PROJECT_COMPLEXITY.out.bam_rmdup \
	| collect //\
	//| DNA_COPYKIT

// USE RNA AND SPATIAL FOR CURIO PIPELINE
	SPATIAL_CURIO(spatial_fq,
					spatial_barcode,
					CELLRANGER_COUNT.out.multiome_outdir)

}

/* See README.md for example run */
