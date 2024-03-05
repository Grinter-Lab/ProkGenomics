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
params.intermedia_files = 'copy'
params.modules = "$baseDir/modules"
params.reference = null
params.software_versions="software_details.txt"
params.adapter_file="TruSeq3-PE.fa"
params.assembly= null
params.version="v.1.0.0" //it has to be one word otherwise will mess up the report
params.report_template="$baseDir/scripts/report.Rmd"
params.logo="$baseDir/scripts/Logo.svg"
params.github="https://github.com/Grinter-Lab/ProkGenomics"
params.default_empty_file="$baseDir/scripts/NO_APPLICABLE"
params.Rrender="$baseDir/scripts/report_render.R"
/*************************************************************************************************************************************************************
 * Channels
 read parameters in the command line 
 ************************************************************************************************************************************************************/

	Channel.value( params.threads )
		.set{ threads_ch }


shortreads = "${params.sample_path}/${params.sample_name}*{1,2}*"
plasmid='plasmids'
chromosome='chromosome'
denovoassembly='denovoassembly'
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
			.map{ file -> tuple (file.baseName,file)}
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
Set defaul channels for omitted steps
*/

myDefaultInputFile = Channel.fromPath(params.default_empty_file)

/*
chInputProcessTwo = chNewInputForProcessTwo.ifEmpty(myDefaultInputFile)

snippy_output = chNewInputForProcessTwo.ifEmpty(myDefaultInputFile)
*/

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
		Type of assembly   : ${params.assembly_type}
		Short reads        : ${shortreads} 
		Long reads         : ${params.longreads}
		Assembly           : ${params.assembly}
		Reference          : ${params.reference}
		Output directory   : ${params.outdir}
		Number of threats  : ${params.threads}
		version            : ${params.version}       
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
include { gtdb } from params.modules

// split assembly
include { split_assembly } from params.modules

// Annotation
include { prokka } from params.modules
include { pharokka } from params.modules

// Comparative Genomics
include { snippy } from params.modules
include { minimap2 } from params.modules

//include { bowtie2 } from params.modules

// Charaterization of genome
//include {prodigal} from params.modules
//include {assembly2feature} from params.modules

//Create a final report of outputs
include { report } from params.modules

/**************************************************************************************************************************************************************
 * workflows subroutine
 **************************************************************************************************************************************************************/


workflow shortreads_QC_workflow{
			fastqc(ch_in_shortreads)
		emit:
			fastqc_html=fastqc.out.html
	}

workflow shortreads_trim_workflow{
		trimmomatic(ch_in_shortreads)
		fastqc(trimmomatic.out.trimmed_reads)
	emit:
		fastqc_trim_html=fastqc.out.html
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
		assemblyqc_path= checkm.out.checkm_tsv
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

workflow split_assembly_workflow{
	take:
		scaffolds
		plasclass_tsv
		checkv_summary
	main:
	 	split_assembly(scaffolds,plasclass_tsv,checkv_summary)
	emit:
		novoassembly_path=split_assembly.out.novoassembly_path
    	chromosome_path=split_assembly.out.chromosome_path
    	plasmid_path=split_assembly.out.plasmid_path.ifEmpty(myDefaultInputFile)
    	phage_path=split_assembly.out.phage_path.ifEmpty(myDefaultInputFile)
}

workflow prokka_scaffolds_workflow{
	take:
		contigs
		denovoassembly
	main:
		prokka(contigs,denovoassembly)
	emit:
		prokka_path=prokka.out.prokka_path
}

workflow prokka_chr_workflow{
	take:
		contigs
		chromosome
	main:
		prokka(contigs,chromosome)
		
	emit:
		prokka_path=prokka.out.prokka_path
}

workflow prokka_plasmids_workflow{
	take:
		contigs
		plasmid
	main:
		prokka(contigs,plasmid)
	emit:
		prokka_path=prokka.out.prokka_path.ifEmpty(myDefaultInputFile)
}

workflow pharokka_workflow{
	take:
		contigs
	main:
		pharokka(contigs)
	emit:
		pharokka_path=pharokka.out.pharokka_path.ifEmpty(myDefaultInputFile)
}




workflow comparative_genomics_workflow{
	take:
	    ch_in_reference
		trimmed_reads
	main:
		snippy(ch_in_reference,trimmed_reads)
		minimap2(ch_in_reference,trimmed_reads)
	emit:
		snippy_path= snippy.out.snippy_path
		minimap2_path= minimap2.out.minimap2_path
}

workflow gtdb_workflow{
	take:
		contigs
	main:
		gtdb(contigs)
	emit:
		gtdb_path= gtdb.out.gtdbtk_summary
}


workflow report_workflow{
	take:
		fastqc_html 
		fastqc_trim_html 
    	novoassembly_path
       	chromosome_path
        plasmid_path
        phage_path
		assembly_qc
        prokka_denovo_path
        prokka_chr_path
        prokka_plasmids_path
        pharokka_path
        plasmid_class
        phage_class
		snippy_path
		gtdb_path
		minimap2_path
       // assembly2gene_table
       // assembly2gene_aligments
       // assembly2gene_peptides
	main:
		report(fastqc_html,
		fastqc_trim_html,
        novoassembly_path,
        chromosome_path,
        plasmid_path,
        phage_path,
		assembly_qc,
        prokka_denovo_path,
        prokka_chr_path,
        prokka_plasmids_path,
        pharokka_path,
        plasmid_class,
        phage_class,
		snippy_path,
		gtdb_path,
		minimap2_path
        //assembly2gene_table,
        //assembly2gene_aligments,
        //assembly2gene_peptides
		)
		
	emit:
		final_report=report.out.report_path
}




/************************************************************************************************************************************************************** 
 * Main workflow
**************************************************************************************************************************************************************/



workflow{

	main:

	if (params.assembly == null){
		if (params.assembly_type=='short'){
			shortreads_QC_workflow()
			shortreads_trim_workflow()
			shortreads_assembly_workflow(shortreads_trim_workflow.out.trimmed_reads)
			if (params.reference )   {comparative_genomics_workflow(ch_in_reference,shortreads_trim_workflow.out.trimmed_reads)
										snippy_output=comparative_genomics_workflow.out.snippy_path
										minimap2_output = comparative_genomics_workflow.out.minimap2_path
										}else{ 
										snippy_output = comparative_genomics_workflow.out.snippy_path.ifEmpty(myDefaultInputFile)
										minimap2_output = comparative_genomics_workflow.out.minimap2_path.ifEmpty(myDefaultInputFile)
										}

		}
		assembly_qc_workflow(shortreads_assembly_workflow.out.scaffolds_path)
		extrachr_workflow(shortreads_assembly_workflow.out.scaffolds)
		split_assembly_workflow(shortreads_assembly_workflow.out.scaffolds,extrachr_workflow.out.plasclass_tsv,extrachr_workflow.out.checkv_summary)
		prokka_chr_workflow(split_assembly_workflow.out.chromosome_path,chromosome)
		prokka_plasmids_workflow(split_assembly_workflow.out.plasmid_path,plasmid)
		prokka_scaffolds_workflow(shortreads_assembly_workflow.out.scaffolds,denovoassembly)
		pharokka_workflow(split_assembly_workflow.out.phage_path)
		gtdb_workflow(split_assembly_workflow.out.chromosome_path)

	}else{
		extrachr_workflow(ch_in_assembly)
		split_assembly_workflow(shortreads_assembly_workflow.out.scaffolds,extrachr_workflow.out.plasclass_tsv,extrachr_workflow.out.checkv_summary)
		prokka_chr_workflow(split_assembly_workflow.out.chromosome_path,chromosome)
		prokka_plasmids_workflow(split_assembly_workflow.out.plasmid_path,plasmid)
		prokka_scaffolds_workflow(shortreads_assembly_workflow.out.scaffolds,denovoassembly)
		pharokka_workflow(split_assembly_workflow.out.phage_path)
		gtdb_workflow(split_assembly_workflow.out.chromosome_path)

	}
	

	fastqc_html_output      = shortreads_QC_workflow.out.fastqc_html.ifEmpty(myDefaultInputFile)
	fastqc_trim_html_output = shortreads_trim_workflow.out.fastqc_trim_html.ifEmpty(myDefaultInputFile)
	scaffolds_output        = shortreads_assembly_workflow.out.scaffolds.ifEmpty(myDefaultInputFile)
	chromosome_path_output  = split_assembly_workflow.out.chromosome_path.ifEmpty(myDefaultInputFile)
	plasmid_path_output     = split_assembly_workflow.out.plasmid_path.ifEmpty(myDefaultInputFile)
	phage_path_output       = split_assembly_workflow.out.phage_path.ifEmpty(myDefaultInputFile)
	prokka_scaffolds_path_output = prokka_scaffolds_workflow.out.prokka_path.ifEmpty(myDefaultInputFile)
	prokka_chr_path         = prokka_chr_workflow.out.prokka_path.ifEmpty(myDefaultInputFile)
	prokka_plasmids_path = prokka_plasmids_workflow.out.prokka_path.ifEmpty(myDefaultInputFile)
	pharokka_path = pharokka_workflow.out.pharokka_path.ifEmpty(myDefaultInputFile)
	plasclass_output = extrachr_workflow.out.plasclass_tsv.ifEmpty(myDefaultInputFile)
    checkv_output = extrachr_workflow.out.checkv_summary.ifEmpty(myDefaultInputFile)
	gtdb_output = gtdb_workflow.out.gtdb_path.ifEmpty(myDefaultInputFile)
	assemblyqc_output = assembly_qc_workflow.out.assemblyqc_path.ifEmpty(myDefaultInputFile)


	report_workflow(	fastqc_html_output,
						fastqc_trim_html_output,
						scaffolds_output,
						chromosome_path_output,
						plasmid_path_output,
						phage_path_output,
						assemblyqc_output,
						prokka_scaffolds_path_output,
						prokka_chr_path,
						prokka_plasmids_path,
						pharokka_path,
						plasclass_output,
        				checkv_output,
        				snippy_output,
						gtdb_output,
						minimap2_output

        				//assembly2gene_table,
        				//assembly2gene_aligments,
        				//assembly2gene_peptides
						)

}



workflow.onComplete { 
	println ( workflow.success ? "\nDone! see the report in ${params.outdir} for more details \n" : "Oops .. something went wrong" )
}
 

