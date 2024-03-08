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
params.default_empty_path="$baseDir/scripts/NO_APPLY/"
params.Rrender="$baseDir/scripts/report_render.R"
/*************************************************************************************************************************************************************
 * Channels
 read parameters in the command line 
 ************************************************************************************************************************************************************/



	Channel.value( params.threads )
		.set{ threads_ch }


	Channel.value( params.sample_name)
		.set{ sample_name}
	

	shortreads = "${params.sample_path}/${params.sample_name}*{1,2}*"

	plasmid        ='plasmids'
	chromosome     ='chromosome'
	denovoassembly ='denovoassembly'



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
Set default channels for omitted steps
*/

/*
touch myDefaultInputFile_assembly
touch myDefaultInputFile_assembly_annotation
touch myDefaultInputFile_chr_annotation
touch myDefaultInputFile_chr_classification
touch myDefaultInputFile_chr_extraction
touch myDefaultInputFile_mapping
touch myDefaultInputFile_phage_annotation
touch myDefaultInputFile_phage_classification
touch myDefaultInputFile_phage_extraction
touch myDefaultInputFile_plasmid_annotation
touch myDefaultInputFile_plasmid_classification
touch myDefaultInputFile_plasmid_extraction
touch myDefaultInputFile_QC_assembly
touch myDefaultInputFile_QC_reads_forward
touch myDefaultInputFile_QC_reads_reverse
touch myDefaultInputFile_QC_trimmed_reads_forward
touch myDefaultInputFile_QC_trimmed_reads_reverse
touch myDefaultInputFile_QC_trimmed_reads_unpaired_forward
touch myDefaultInputFile_QC_trimmed_reads_unpaired_reverse
touch myDefaultInputFile_SNV_detection
*/


Channel
 	.fromPath(params.default_empty_path)
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile}

//myDefaultInputFile.view()



Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_assembly")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_assembly}

//myDefaultInputFile_assembly.view()


Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_assembly_annotation")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_assembly_annotation}

//myDefaultInputFile_assembly_annotation.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_chr_annotation")
	.map { file -> tuple (file.baseName,file)}
	.set{myDefaultInputFile_chr_annotation}

//myDefaultInputFile_chr_annotation.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_chr_classification")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_chr_classification}

//myDefaultInputFile_chr_classification.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_chr_extraction")
	.map { file -> tuple (file.baseName,file)}
	.set{myDefaultInputFile_chr_extraction}

//myDefaultInputFile_chr_extraction.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_mapping")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_mapping}

//myDefaultInputFile_mapping.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_phage_annotation")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_phage_annotation}

//myDefaultInputFile_phage_annotation.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_phage_classification")
	.map { file -> tuple (file.baseName,file)}
	.set{myDefaultInputFile_phage_classification}

//myDefaultInputFile_phage_classification.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_phage_extraction")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_phage_extraction}

//myDefaultInputFile_phage_extraction.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_plasmid_annotation")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_plasmid_annotation}

//myDefaultInputFile_plasmid_annotation.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_plasmid_classification")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_plasmid_classification}

//myDefaultInputFile_plasmid_classification.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_plasmid_extraction")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_plasmid_extraction}

//myDefaultInputFile_plasmid_extraction.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_QC_assembly")
	.map { file -> tuple (file.baseName,file)}
	.set{myDefaultInputFile_QC_assembly}

//myDefaultInputFile_QC_assembly.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_QC_reads*")
	.collect()
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_QC_reads}

//myDefaultInputFile_QC_reads.view()

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_QC_trimmed_*")
	.collect()
	.map { file -> tuple (file.baseName,file)}
	.set{myDefaultInputFile_QC_trimmed_reads}

//myDefaultInputFile_QC_trimmed_reads.view()


Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_SNV_detection")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_SNV_detection}



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
    	plasmid_path=split_assembly.out.plasmid_path
    	phage_path=split_assembly.out.phage_path
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
		prokka_path=prokka.out.prokka_path
}

workflow pharokka_workflow{
	take:
		contigs
	main:
		pharokka(contigs)
	emit:
		pharokka_path=pharokka.out.pharokka_path
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
		sample_name
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
	 report(
		sample_name,
		fastqc_html,
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
				if (params.reference ){ comparative_genomics_workflow(ch_in_reference,shortreads_trim_workflow.out.trimmed_reads)
										snippy_output = comparative_genomics_workflow.out.snippy_path
										minimap2_output = comparative_genomics_workflow.out.minimap2_path
										}else{ 
										snippy_output = myDefaultInputFile_SNV_detection
										minimap2_output = myDefaultInputFile_mapping
										}

		
			assembly_qc_workflow(shortreads_assembly_workflow.out.scaffolds_path)
			extrachr_workflow(shortreads_assembly_workflow.out.scaffolds)
			split_assembly_workflow(shortreads_assembly_workflow.out.scaffolds,extrachr_workflow.out.plasclass_tsv,extrachr_workflow.out.checkv_summary)
			prokka_chr_workflow(split_assembly_workflow.out.chromosome_path,chromosome)
			prokka_plasmids_workflow(split_assembly_workflow.out.plasmid_path,plasmid)
			prokka_scaffolds_workflow(shortreads_assembly_workflow.out.scaffolds,denovoassembly)
			pharokka_workflow(split_assembly_workflow.out.phage_path)
			gtdb_workflow(split_assembly_workflow.out.chromosome_path)


			fastqc_html_output = shortreads_QC_workflow.out.fastqc_html.ifEmpty{ myDefaultInputFile_QC_reads }
			fastqc_trim_html_output = shortreads_trim_workflow.out.fastqc_trim_html.ifEmpty{myDefaultInputFile_QC_trimmed_reads }
			scaffolds_output = shortreads_assembly_workflow.out.scaffolds.ifEmpty{myDefaultInputFile_assembly }
			assemblyqc_output = assembly_qc_workflow.out.assemblyqc_path.ifEmpty{myDefaultInputFile_QC_assembly }
			chromosome_path_output = split_assembly_workflow.out.chromosome_path.ifEmpty{myDefaultInputFile_chr_extraction }
			plasmid_path_output = split_assembly_workflow.out.plasmid_path.ifEmpty{myDefaultInputFile_plasmid_extraction }
			phage_path_output = split_assembly_workflow.out.phage_path.ifEmpty{myDefaultInputFile_phage_extraction }
			prokka_scaffolds_path_output = prokka_scaffolds_workflow.out.prokka_path.ifEmpty{myDefaultInputFile_assembly_annotation }
			prokka_chr_path = prokka_chr_workflow.out.prokka_path.ifEmpty{myDefaultInputFile_chr_annotation }
			prokka_plasmids_path = prokka_plasmids_workflow.out.prokka_path.ifEmpty{myDefaultInputFile_plasmid_annotation }
			pharokka_path = pharokka_workflow.out.pharokka_path.ifEmpty{myDefaultInputFile_phage_annotation }
			plasclass_output = extrachr_workflow.out.plasclass_tsv.ifEmpty{myDefaultInputFile_plasmid_classification }
			checkv_output = extrachr_workflow.out.checkv_summary.ifEmpty{myDefaultInputFile_phage_classification }
			gtdb_output = gtdb_workflow.out.gtdb_path.ifEmpty{myDefaultInputFile_chr_classification }

			}else if(params.assembly_type=='long')
				{

			}else if(params.assembly_type=='hybrid')
				{

					}

		}else{

			fastqc_html_output = myDefaultInputFile_QC_reads
			fastqc_trim_html_output = myDefaultInputFile_QC_trimmed_reads
			scaffolds_output = myDefaultInputFile_assembly
			assemblyqc_output = myDefaultInputFile_QC_assembly
			snippy_output = myDefaultInputFile_SNV_detection
			minimap2_output = myDefaultInputFile_mapping

			extrachr_workflow(ch_in_assembly)
			split_assembly_workflow(ch_in_assembly,extrachr_workflow.out.plasclass_tsv,extrachr_workflow.out.checkv_summary)
			prokka_chr_workflow(split_assembly_workflow.out.chromosome_path,chromosome)
			prokka_plasmids_workflow(split_assembly_workflow.out.plasmid_path,plasmid)
			prokka_scaffolds_workflow(ch_in_assembly,denovoassembly)
			pharokka_workflow(split_assembly_workflow.out.phage_path)
			gtdb_workflow(split_assembly_workflow.out.chromosome_path)

			chromosome_path_output = split_assembly_workflow.out.chromosome_path.ifEmpty{ myDefaultInputFile_chr_extraction }
			plasmid_path_output = split_assembly_workflow.out.plasmid_path.ifEmpty{myDefaultInputFile_plasmid_extraction }
			phage_path_output = split_assembly_workflow.out.phage_path.ifEmpty{myDefaultInputFile_phage_extraction }
			prokka_scaffolds_path_output = prokka_scaffolds_workflow.out.prokka_path.ifEmpty{myDefaultInputFile_assembly_annotation }
			prokka_chr_path = prokka_chr_workflow.out.prokka_path.ifEmpty{myDefaultInputFile_chr_annotation }
			prokka_plasmids_path = prokka_plasmids_workflow.out.prokka_path.ifEmpty{myDefaultInputFile_plasmid_annotation }
			pharokka_path = pharokka_workflow.out.pharokka_path.ifEmpty{myDefaultInputFile_phage_annotation }
			plasclass_output = extrachr_workflow.out.plasclass_tsv.ifEmpty{myDefaultInputFile_plasmid_annotation }
			checkv_output = extrachr_workflow.out.checkv_summary.ifEmpty{myDefaultInputFile_phage_classification }
			gtdb_output = gtdb_workflow.out.gtdb_path.ifEmpty{myDefaultInputFile_chr_classification }
		}



	report_workflow( sample_name,
					 fastqc_html_output,
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
					 minimap2_output,
					 gtdb_output

        				//assembly2gene_table,
        				//assembly2gene_aligments,
        				//assembly2gene_peptides
					)
}




workflow.onComplete { 
	println ( workflow.success ? "\nDone! see the report in ${params.outdir} for more details \n" : "Oops .. something went wrong" )
}
 

