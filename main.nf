#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/**************************************************************************************************************************************************************

# Running at local computer: 
cd /Users/lper0012/Dropbox/BioPlatform/NextFlow_ProkAssembly/pipeline
export JAVA_HOME='/opt/homebrew/opt/openjdk@11/'
nextflow run main.nf --assembly_type 'short' --shortreads "$baseDir/test/*{R1,R2}.fastq"  --outdir "$baseDir/../results/" --modules "$baseDir/modules"
nextflow run main.nf --assembly_type 'short' --sample_name '1-77321' -profile conda

# Running at bio1  
cd /home/lper0012/tasks/rhys.grinter/pipeline
nextflow run main.nf --assembly_type 'short' --sample_name '1-77321' -profile conda


* Out line of the pipeline
 Pipeline 1- Assembly and annotation pipeline:
	QC reads
	trimming
	Assembly - de novo (unicycler) -either hybrid, short or long
	QC assembly quality (Busco)
	Identification of the species (GTDB)
	plasmids (plasclass)
	Phage prediction/ prophage (DBSCAN-SWA) PhiSpy 
	Phage matching CheckV
	Annotation (prokka for prokaryotes or pharokka for phage)
	Assembly to ref - coverage, core genome, variants, mutations,SNPs, Rearrangements and larger deletions(Diversitools/Parsnip/clinker)
		-Output as fasta, gff, genbank files and various summary tables 

Pipeline 2: Assembly to gene pipeline: 
	Create a local blast database using provided proteins (ie BamA) for all new annotations to be compared against.
other:
		Produce documentation for using pipelines and give a tutorial to the team on how to use the pipelines.
*
**************************************************************************************************************************************************************/

/**************************************************************************************************************************************************************
 * pipeline input parameters
 **************************************************************************************************************************************************************/

// make people input folder input and prefix of the files only and create shortreads with that
params.sample_path="$PWD/rawdata/"
params.sample_name=''
params.assembly_type = 'short'
params.longreads = null
params.threads = (16)
params.outdir = "$PWD/${params.sample_name}_results"
params.publish_dir_mode = 'copy'
params.modules = "$baseDir/modules"
params.reference = null
params.software_versions="software_details.txt"
params.adapter_file="TruSeq3-PE.fa"
params.assembly=null

/*************************************************************************************************************************************************************
 * Channels
 read parameters in the command line 
 ************************************************************************************************************************************************************/

	Channel.value( params.threads )
		.set{ threads_ch }


shortreads = "${params.sample_path}/${params.sample_name}*{1,2}*"
/*
Pick the type of de novo-assembly process
*/
	if (params.assembly_type == 'hybrid') {
		Channel.fromFilePairs(shortreads, checkIfExists: true)
				.set{  ch_in_shortreads }  
		Channel.fromPath(params.longreads, checkIfExists: true)
				.map{ file -> tuple (file.baseName,file)}
				.set{ch_in_longreads}
		}else if (params.assembly_type == 'short') {
				Channel.fromFilePairs(shortreads, checkIfExists: true)
					.set{ ch_in_shortreads}
		}else if (params.assembly_type == 'long') {    
				Channel.fromPath(params.longreads, checkIfExists: true)
				.map{ file -> tuple (file.baseName,file)}
				.set{ch_in_longreads}
		}
	
/*
If reference provided read the parameter. If reference not provided then find a good reference in next steps
*/

	if (params.reference){      
		Channel.fromPath(params.reference, checkIfExists: true)
			.set{ch_in_reference}
	}

/*
If the assembly has been done provide assembly file
*/

	if (params.assembly){      
		Channel.fromPath(params.assembly, checkIfExists: true)
			.map{ file -> tuple (file.baseName,file)}
			.set{ch_in_assembly}		
	}

 	
/*
Define a env or container
*/
 	//if(params.enable_conda){ch_in_env="conda"}
	
/**************************************************************************************************************************************************************
 * Print parameters
 **************************************************************************************************************************************************************/


log.info """\
		ProkGenomics: Prokaryotic Genomics Pipeline
		Author: Laura Perlaza-Jimenez (PhD)
		Rhys Grinter Laboratory
		======================================================================
		Type of Assembly   : ${params.assembly_type}
		Short Reads        : ${shortreads} 
		Long Reads         : ${params.longreads}
		Assembly           : ${params.assembly}
		Reference          : ${params.reference}
		Output directory   : ${params.outdir}
		Number of Threats  : ${params.threads}
		======================================================================        


         """
         .stripIndent()



/**************************************************************************************************************************************************************
* Include the following modules from modules file
**************************************************************************************************************************************************************/

// Reads QC
include { fastqc } from params.modules
include { trimmomatic } from params.modules
//include { fastxtoolkit } from params.modules

// De-novo assembly
include { unicycler_short } from params.modules

// Assembly QC
//include { busco } from params.modules
include { checkm } from params.modules

// Assembly classification of contigs
include { checkv } from params.modules
include { plasclass } from params.modules

// Taxonomy classification
//include { GTDB } from params.modules

// Annotation
include { prokka } from params.modules
include { pharokka } from params.modules

// Comparative Genomics
include { snippy } from params.modules
//include { minimap2 } from params.modules

// Charaterization of genome
//include {assembly2feature} from params.modules

/**************************************************************************************************************************************************************
 * workflows subroutine
 **************************************************************************************************************************************************************/


workflow shortreads_QC_workflow{
		fastqc(ch_in_shortreads)

		//emit:
		//assembly = metaspades.out.assemblyshort
		//ch_in_reads=ch_in_shortreads  
	}

workflow shortreads_trim_workflow{
		trimmomatic(ch_in_shortreads)
		fastqc(trimmomatic.out.trimmed_reads)
	emit:
		trimmed_reads=trimmomatic.out.trimmed_reads
}

workflow shortreads_assembly_workflow{
	take:
		trimmed_reads
	main:
		unicycler_short(trimmed_reads)
	emit:
		scaffolds_path=unicycler_short.out.scaffolds_path
		scaffolds=unicycler_short.out.scaffolds
}

workflow assembly_qc_workflow{
	take:
		scaffolds_path
	main:
		checkm(scaffolds_path)
	emit:
		assemblyQC=checkm.out.checkm_tsv
}


workflow extrachr_workflow{
	take:
		scaffolds
	main:
		plasclass(scaffolds)
		checkv(scaffolds)
	emit:
		plasclass_tsv=plasclass.out.plasclass_tsv
		checkv_summary=checkv.out.checkv_summary
}

workflow annotation_workflow{
	take:
		scaffolds
	main:
		prokka(scaffolds)
		pharokka(scaffolds)

}

workflow SNV_workflow{
		shortreads=ch_in_shortreads
		referencegenome=ch_in_reference
	main:
		snippy(referencegenome,shortreads)
}

/************************************************************************************************************************************************************** 
 * Main workflow
**************************************************************************************************************************************************************/



workflow{

	main:

	if (params.assembly == null){
		if (params.assembly_type=='short'){
			if ( params.reference ) { SNV_workflow()}
			shortreads_QC_workflow()
			shortreads_trim_workflow()
			shortreads_assembly_workflow(shortreads_trim_workflow.out.trimmed_reads)
		}

		if (params.assembly_type=='long'){
			if ( params.reference ) { SNV_workflow()}
			shortreads_QC_workflow()
			shortreads_trim_workflow()
			shortreads_assembly_workflow(shortreads_trim_workflow.out.trimmed_reads)
		}
		assembly_qc_workflow(shortreads_assembly_workflow.out.scaffolds_path)
		extrachr_workflow(shortreads_assembly_workflow.out.scaffolds)
		annotation_workflow(shortreads_assembly_workflow.out.scaffolds)
	}
	//scaffolds_path=ch_in_assembly
	//assembly_qc_workflow(ch_in_assembly_path)
	extrachr_workflow(ch_in_assembly)
	//annotation_workflow(ch_in_assembly)
		



//	assembly.view()
//	mapping_workflow(assembly,ch_in_reads)
//	mappedsam=mapping_workflow.out[0]
//	mappedbam=mapping_workflow.out[1]
//	binning_workflow(assembly,mappedsam,mappedbam)
//	refbinning=binning_workflow.out[0]
//	dRep_GTDB_workflow(refbinning)

}


workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir\n" : "Oops .. something went wrong" )
}
 


































