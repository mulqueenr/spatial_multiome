//Nextflow pipeline for processing Navin lab spatial multiome//
// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
// DNA
params.dna_flowcellDir = "/volumes/seq/flowcells/MDA/nextseq2000/2025/250227_RM_CurioWGS_scalemet" //Sequencing run flowcell dir
params.dna_samplesheet = "DNA_SampleSheet.csv"
params.spatial_barcode = "/volumes/USR3/ctang4/10X_processing_RNA/Ryan_Trekker/022525_CuioWGS_DCIS41T/U0028_003_BeadBarcodes.txt"

// RNA
params.rna_flowcellDir = "/Volumes/seq/flowcells/MDA/nextseq2000/2025/250220_RM_CuioWGS_RNA" //Sequencing run flowcell dir
params.rna_samplesheet = "RNA_SimpleSampleSheet.csv"

//REF
params.ref="/volumes/USR2/Ryan/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
params.cellranger_arc="/volumes/USR2/Ryan/tools/cellranger-arc-2.0.2/cellranger-arc"
params.cellranger_atac="/volumes/USR2/Ryan/tools/cellranger-atac-2.1.0/cellranger-atac"
params.cellranger_rna="/volumes/USR2/Ryan/tools/cellranger-9.0.1/cellranger"
params.max_cpus="50"

//output
params.outname = "250129_spatialdna"
params.outdir = "/volumes/USR2/Ryan/projects/spatial_wgs/data/250129_First_Experiment2"
params.date = "20250226"
//library parameters
params.cell_try="5000" //Based on expected cell count from library generation
params.samplesheet="/volumes/USR2/Ryan/projects/spatial_wgs/data/250129_First_Experiment/DNA_SampleSheet.csv" //Based on expected cell count from library generation

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
		Cellranger ARC install : ${params.cellranger}
		Max cpus : ${params.max_cpus}
		================================================

""".stripIndent()

// BCL TO FASTQ PIPELINE FOR GENERATING SINGLE-CELL FASTQs
process DNA_BCL_TO_FASTQ { 
	//Generate Undetermined Fastq Files from BCL Files.
    //Count GEM indexes and generate a white list for splitting
	//Assumes Y151;I10;U16;Y151 sequencing cycles unless specified as input parameter
	//bcl-convert requires write access to "/var/logs/bcl-convert", so we just bind a dummy one
	cpus "${params.max_cpus}"
	input:
		path(dna_flowcellDir)
		path(dna_samplesheet)
	output:
		path("*fastq.gz")
    script:
		"""
        #Run initial bcl convert and count gem indexes to determine whitelist for splitting
        task_cpus=\$(expr ${task.cpus} / 3)

        bcl-convert \\
        --bcl-input-directory ${dna_flowcellDir} \\
        --bcl-num-conversion-threads \$task_cpus \\
        --bcl-num-compression-threads \$task_cpus \\
        --bcl-num-decompression-threads \$task_cpus \\
		--bcl-only-matched-reads true \\
        --sample-sheet ${dna_samplesheet} \\
        --no-lane-splitting true \\
        --output-directory . \\
        --force
		"""
}

process RNA_CELLRANGER_MKFASTQ{
	//Run cellranger on RNA samples, this works straight from bcl files.
	// Run mkfastq then run count
	cpus "${params.max_cpus}"

	input:
		path(rna_flowcellDir)
		path(rna_samplesheet)
	output:
		path("${params.outname}/outs/fastq_path/*/${params.outname}_rna*{I1,I2,R1,R2}_001.fastq.gz"), emit: transcriptome
		path("${params.outname}/outs/fastq_path/*/${params.outname}_spatial*{I1,I2,R1,R2}_001.fastq.gz"), emit: spatial

    script:
		"""
		${params.cellranger_rna} mkfastq \\
		--run=${rna_flowcellDir} \\
		--id=${params.outname} \\
		--samplesheet=${rna_samplesheet} \\
		--localcores=${params.max_cpus} \\
		--delete-undetermined \\
		--localmem=300
		"""

}

process RNA_CELLRANGER_COUNT {
	//Run cellranger on DNA samples, to generate GEM-indexed bam file.
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/rna_cellranger", mode: 'copy', overwrite: true, pattern: "./outs/*"

	input:
		path(rna_fqDir), stageAs: 'rna_fq/*'

	output:
		path("./outs/possorted_bam.bam"), emit: bam
		path("./outs/*"), emit: outdir

    script:
		"""
      	${params.cellranger_rna} count \\
		--fastqs="\${PWD}/rna_fq/" \\
		--reference=${params.ref} \\
		--id=${params.outname} \\
		--chemistry=ARC-v1 \\
		--localcores=${params.max_cpus} \\
        --localmem=300
		"""
}

process DNA_CELLRANGER_COUNT {
	//Run cellranger on DNA samples, to generate GEM-indexed bam file.
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/dna_cellranger", mode: 'copy', overwrite: true, pattern: "./outs/*"

	input:
		path(dna_fqDir), stageAs: 'dna_fq/*'

	output:
		path("./outs/possorted_bam.bam"), emit: dna_bam
		path("./outs/*"), emit: dna_outdir

    script:
		"""
        ${params.cellranger_atac} count \\
		--fastqs="\${PWD}/dna_fq/" \\
		--reference=${params.ref} \\
		--id=${params.outname} \\
		--sample=${params.outname}_dna \\
		--chemistry=ARC-v1 \\
		--localcores=${params.max_cpus} \\
        --localmem=300
		"""
}
process SPATIAL_CURIO {
	//Run Curio Trekker Pipeline to generate spatial location
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/spatial", mode: 'copy', overwrite: true

	input:
		tuple path(fq_i1),path(fq_i2),path(fq_r1),path(fq_r2)
		path(spatial_barcode)

	output:
		path("./outs/possorted_bam.bam"), emit: bam
		path("./outs/*"), emit: outdir

    script:
		"""
		sample_name="${params.outname}_rna"
		experiment_date="${params.date}"

		echo 'sample,sc_sample,experiment_date,barcode_file,fastq_1,fastq_2,sc_outdir,sc_platform,profile,subsample,cores' > samplesheet.trekker.csv
		echo "\${sample_name},\${sample_name},\${experiment_dat},${spatial_barcode},${fq_r1},${fq_r2},\${PWD},TrekkerU_C,singularity,no,${task.cpus}" >> samplesheet.trekker.csv

		bash /volumes/USR2/Ryan/tools/curiotrekker-v1.1.0/nuclei_locater_toplevel.sh \\
		samplesheet.trekker.csv
		"""
}


process DNA_SPLIT_BAM {
	//Use timoast sinto to split
	//conda install sinto

	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/dna_cellranger/sc_bam", mode: 'copy', overwrite: true 

	input:
		tuple path(dna_bam),path(dna_barcodes)

	output:
		path("./outs/possorted_bam.bam"), emit: bam
		path("./outs/*"), emit: outdir

    script:
		"""
		sinto filterbarcodes \\
		--bam ${bam} \\
		--cells \\
		-p ${task.cpus} \\
		--barcodetag "CB" \\
		--outdir ./sc_bam
		"""
}


workflow {
// BCL TO FASTQ PIPELINE FOR SPLITTING FASTQS		
	dna_flowcell_dir = Channel.fromPath(params.dna_flowcellDir)
	dna_samplesheet = Channel.fromPath(params.dna_samplesheet)
	rna_flowcell_dir = Channel.fromPath(params.rna_flowcellDir)
	rna_samplesheet = Channel.fromPath(params.rna_samplesheet)
	spatial_barcode = Channel.fromPath(params.spatial_barcode)

	//Generate DNA Fastqs
	dna_out = \
	DNA_BCL_TO_FASTQ(dna_flowcell_dir,dna_samplesheet) \
	| collect \
	| DNA_CELLRANGER_COUNT

	//Generate RNA Fastqs and split
	rna_fq_in = \
	RNA_CELLRANGER_MKFASTQ(rna_flowcell_dir,rna_samplesheet)

	rna_out = rna_fq_in.transcriptome \
	| collect \
	| RNA_CELLRANGER_COUNT

	//Run cellranger on DNA/RNA

	//Run curio pipeline on spatial
	SPATIAL_CURIO(rna_fq_in.spatial,spatial_barcode)
	//Output DNA BAM goes into copykit
	//Output RNA matrix goes into seurat object
	//Output RNA Spatial goes into metadata for seurat object

}

/* See README.md for example run */
