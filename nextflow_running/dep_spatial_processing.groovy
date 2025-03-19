//Nextflow pipeline for processing Navin lab spatial multiome//
// Declare syntax version
nextflow.enable.dsl=2

// Script parameters
// DNA
params.dna_flowcellDir = "/home/rmulqueen/projects/spatial_wgs/seq/250227_RM_CurioWGS_scalemet" //Sequencing run flowcell dir
params.dna_samplesheet = "DNA_SampleSheet.csv"

params.spatial_barcode = "/home/rmulqueen/tools/curiotrekker-v1.1.0/U0028_003_BeadBarcodes.txt"

// RNA
params.rna_flowcellDir = "/home/rmulqueen/projects/spatial_wgs/seq/250220_RM_CuioWGS_RNA" //Sequencing run flowcell dir
params.rna_samplesheet = "RNA_SimpleSampleSheet.csv"

//REF
params.ref="/home/rmulqueen/ref/refdata-cellranger-arc-GRCh38-2020-A-2.0.0"
params.src="/home/rmulqueen/projects/spatial_wgs/tools/spatial_multiome/src"
params.cellranger_arc="/home/rmulqueen/tools/cellranger-arc-2.0.2/cellranger-arc"
params.cellranger_atac="/home/rmulqueen/tools/cellranger-atac-2.1.0/cellranger-atac"
params.cellranger_rna="/home/rmulqueen/tools/cellranger-9.0.1/cellranger"
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
		Cellranger ATAC install : ${params.cellranger_atac}
		Cellranger RNA install : ${params.cellranger_rna}

		Max cpus : ${params.max_cpus}
		================================================

""".stripIndent()

// BCL TO FASTQ PIPELINE FOR GENERATING SINGLE-CELL FASTQs
process DNA_BCL_TO_FASTQ { 
	//Generate Undetermined Fastq Files from BCL Files.
	//bcl-convert requires write access to "/var/logs/bcl-convert", so we just bind a dummy one if we add a sif
	cpus "${params.max_cpus}"
	containerOptions "--bind ${params.outdir}/logs:/var/log/bcl-convert"	
	label 'amethyst'

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

process DNA_CELLRANGER_COUNT {
	//Run cellranger on DNA samples, to generate GEM-indexed bam file.
	//This throws an error when done, not sure why
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/dna_cellranger", mode: 'copy', overwrite: true

	input:
		path(dna_fqDir), stageAs: 'dna_fq/*'

	output:
		tuple path("./${params.outname}/outs/possorted_bam.bam"), path("./${params.outname}/outs/possorted_bam.bam.bai"),path("./${params.outname}/outs/filtered_peak_bc_matrix/barcodes.tsv"), emit: dna_bam
		path("./${params.outname}/outs/*"), emit: dna_outdir

    script:
		"""
        ${params.cellranger_atac} count \\
		--fastqs="\${PWD}/dna_fq/" \\
		--reference=${params.ref} \\
		--id=${params.outname} \\
		--localcores=${params.max_cpus} \\
		--chemistry=ARC-v1 \\
        --localmem=1000
		"""

}

process DNA_SPLIT_BAM {
	//Use timoast sinto to split bam by CB tag
	//conda install sinto
	//cpu limit hard set because this causes a lot of i/o
	cpus 50
	publishDir "${params.outdir}/dna_cellranger/sc_dna_bam", mode: 'copy', overwrite: true 

	input:
		tuple path(dna_bam),path(dna_bam_bai),path(dna_barcodes)

	output:
		path("./sc_dna_bam/*bam"), emit: bam

    script:
		"""
		#make splitting barcode list from whitelist
		awk -F, 'OFS="\t" {print \$1,\$1}' ${dna_barcodes} > cell_id.tsv

		#split to chunks of 500 cells for i/o purposes
		split -l 500 --numeric-suffixes cell_id.tsv cell_id.split.

		#run cell splitting for each 500 chunk
		for i in cell_id.split* ; do
		sinto filterbarcodes --bam ${bam} --cells \$i -p ${task.cpus} --barcodetag "CB" --outdir ./sc_dna_bam ;
		done
		"""
}


process DNA_PROJECT_COMPLEXITY {
	//Use picard tools to project library complexity
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/dna_cellranger/sc_bam_dedup", mode: 'copy', overwrite: true, pattern: "*rmdup.bam"
	publishDir "${params.outdir}/reports/dna/complexity", mode: 'copy', overwrite: true , pattern: "*metrics.txt"
	publishDir "${params.outdir}/reports/dna/rmdup", mode: 'copy', overwrite: true , pattern: "*rmdup.stats.txt"

	input:
		path(bam)

	output:
		path("*rmdup.bam"), emit: bam_rmdup
		path("*complex_metrics.txt"), emit: complexity_metrics
		path("*rmdup.stats.txt"), emit: rmdup_metrics

    script:
	/*function proj_complexity() {
	bam=$1
	bam_simplename=${bam::-4}
	samtools sort -T . -n -o - ${bam} | \
	samtools fixmate -m - - | \
	samtools sort -T . -o - - | \
	samtools markdup -s - ${bam_simplename}.rmdup.bam 2> ${bam_simplename}.rmdup.stats.txt

	java -jar ~/tools/picard.jar \
	EstimateLibraryComplexity \
	MAX_OPTICAL_DUPLICATE_SET_SIZE=-1 \
	I=${bam_simplename}.rmdup.bam \
	O=${bam_simplename}.complex_metrics.txt
	}

	export -f proj_complexity
	parallel -j 100 proj_complexity ::: $(ls *-1.bam)

	for i in *complex_metrics.txt;
	do grep "^Unknown" $i | awk -v cellid=${i::-20} 'OFS="," {print cellid,$3,$9,$10}' ; done > dna_cell_metrics.csv
	*/
		"""
		"""
}

process DNA_COPYKIT {
	// CNV PROFILING 
	//COPYKIT FOR CLONE CALLING BY CNVS
	//Run CopyKit and output list of bam files by clones
	cpus "${params.max_cpus}"
	label 'cnv'
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"
	publishDir "${params.outdir}/cnv_calling", mode: 'copy', pattern: "*{tsv,rds}"
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
	publishDir "${params.outdir}/rna_cellranger", mode: 'copy', overwrite: true

	input:
		path(rna_fqDir), stageAs: 'rna_fq/*'

	output:
		path("*"), emit: outdir

    script:
		"""
      	${params.cellranger_rna} count \\
		--fastqs="\${PWD}/rna_fq/" \\
		--transcriptome=${params.ref} \\
		--id=${params.outname} \\
		--create-bam=true \\
		--chemistry=ARC-v1 \\
		--localcores=${params.max_cpus} \\
		--localmem=300 2> cellranger.log
		"""
}

process RNA_SEURAT_OBJECT_GENERATION {
	//Make Seurat object from cellranger output
	cpus "${params.max_cpus}"
	label 'amethyst' //I havent made a dedicated sif for spatial project yet.
	publishDir "${params.outdir}/plots", mode: 'copy', overwrite: true, pattern: '*pdf'
	containerOptions "--bind ${params.src}:/src/,${params.outdir}"

	input:
		path(outdir)

	output:
		path("*.seuratObj.rds"), emit: rds
		path("*pdf"), emit: seurat_plots


    script:
		"""
      	/src/seurat_cellranger_output.R \\
		--input_dir ./${params.outname} \\
		--out_name ${params.outname}
		"""
	
}

process SPATIAL_CURIO {
	//Run Curio Trekker Pipeline to generate spatial location
	cpus "${params.max_cpus}"
	publishDir "${params.outdir}/spatial", mode: 'copy', overwrite: true

	input:
		tuple path(fq_i1),path(fq_i2),path(fq_r1),path(fq_r2)
		path(spatial_barcode)
		path(sc_outdir), stageAs: 'sc_outdir/*'

	output:
		path("*")

    script:
	"""
	sample_name="${params.outname}_rna"
	cp \$(realpath ./sc_outdir/${params.outname}/outs/filtered_feature_bc_matrix) ./filtered_feature_bc_matrix 

	echo 'sample,sc_sample,experiment_date,barcode_file,fastq_1,fastq_2,sc_outdir,sc_platform,profile,subsample,cores' > samplesheet.trekker.csv
	echo "${sample_name},${sample_name},${params.date},${spatial_barcode},${fq_r1},${fq_r2},\${PWD}/filtered_feature_bc_matrix,TrekkerU_C,singularity,no,${task_cpus}" >> samplesheet.trekker.csv

	bash /home/rmulqueen/tools/curiotrekker-v1.1.0/nuclei_locater_toplevel.sh \\
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

//Generate copy number calls from DNA data
	DNA_BCL_TO_FASTQ(dna_flowcell_dir,dna_samplesheet) //
	| collect \
	| DNA_CELLRANGER_COUNT

	
	DNA_CELLRANGER_COUNT.out.dna_bam \
	| DNA_SPLIT_BAM \
	| DNA_PROJECT_COMPLEXITY

	DNA_PROJECT_COMPLEXITY.out.bam_rmdup \
	| collect \
	| DNA_COPYKIT
	

	//Generate seurat object from RNA data
	RNA_CELLRANGER_MKFASTQ(rna_flowcell_dir,rna_samplesheet)

	RNA_CELLRANGER_MKFASTQ.out.transcriptome \
	| collect \
	| RNA_CELLRANGER_COUNT

	//RNA_SEURAT_OBJECT_GENERATION(RNA_CELLRANGER_COUNT.out.outdir)

	//Generate spatial information from curio oligoes
	SPATIAL_CURIO(RNA_CELLRANGER_MKFASTQ.out.spatial,spatial_barcode,RNA_CELLRANGER_COUNT.out.outdir)
}

/* See README.md for example run */
