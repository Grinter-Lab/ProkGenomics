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
nextflow run main.nf --assembly_type 'short' --sample_name '1-77321' -profile singularity --reference unicycler/1-77321-LFA6246_L2/assembly.fasta -resume


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
params.sample_name=null
params.assembly_type = 'short'
params.longreads = null
params.threads = (16)
params.outdir = "$PWD/${params.sample_name}_results"
params.publish_dir_mode = 'copy'
params.publish_intermediate_files  = 'copy'
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
params.run_classification="TRUE"
params.db_gtdb_path="$baseDir//DB/db_gtdb/release214/"
params.genes_interest=null
params.percentage=80
params.keep_intermediate_files=false
params.help=false
params.cleanup=true

/*************************************************************************************************************************************************************
 *Print
 ************************************************************************************************************************************************************/


if( params.help | params.sample_name == null) {

log.info """
·········································································································
:  ██████╗ ██████╗  ██████╗ ██╗  ██╗ ██████╗ ███████╗███╗   ██╗ ██████╗ ███╗   ███╗██╗ ██████╗███████╗  :
:  ██╔══██╗██╔══██╗██╔═══██╗██║ ██╔╝██╔════╝ ██╔════╝████╗  ██║██╔═══██╗████╗ ████║██║██╔════╝██╔════╝  :
:  ██████╔╝██████╔╝██║   ██║█████╔╝ ██║  ███╗█████╗  ██╔██╗ ██║██║   ██║██╔████╔██║██║██║     ███████╗  :
:  ██╔═══╝ ██╔══██╗██║   ██║██╔═██╗ ██║   ██║██╔══╝  ██║╚██╗██║██║   ██║██║╚██╔╝██║██║██║     ╚════██║  :
:  ██║     ██║  ██║╚██████╔╝██║  ██╗╚██████╔╝███████╗██║ ╚████║╚██████╔╝██║ ╚═╝ ██║██║╚██████╗███████║  :
:  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝ ╚═╝     ╚═╝╚═╝ ╚═════╝╚══════╝  :
·········································································································
ProkGenomics: Prokaryotic Genomics Pipeline
Author: Laura Perlaza-Jimenez (PhD) ~ Monash Genomics and Bioinformatics Platform
Rhys Grinter Laboratory
======================================================================================================================================================
ERROR: Please provide the minimum options 
=============================================
Usage:
    ProkGenomics --sample_name <sample_name> -profile singularity
	
Input:
--sample_path           	The default path for the reads is the folder rawdata in the working directory. if you have your reads somewhere 
							else you should set this parameter to that path.

--sample_name           	The sample name is the prefix of your samples files. it doesn't have a default because I don't know your sample 
							names. Please don't use sample names with spaces in them. Best approach is to use the name of the file as it comes 
							from the sequencing facility

--assembly_type         	This parameter can be short long or hybrid. The default is 'short'. if you have short reads you don't have to 
							specify this parameter. If you pick the argument long or hybrid the longreads parameter should be specify. For 
							hybrid make sure to give a path for short and long reads. 

--longreads     			Path to the long reads. 

--threads               	Number of threats to use. More threats faster your processing. Make sure you know what is available for you.

--outdir                	The results will be in a folder in the working directory with the same sample name and _results ex. 1-77321_results.

--assembly					If you provide a path to an assembly from the reads, the assembly steps will be skipped.

--reference 				If you have a reference genome put the path here. This will activate all the comparative genomics steps. This file 
							can be formatted as FASTA or GENBANK. If you provide a GENBANK file your Single Nucleotide Variant file will be annotated 
							(tell you what gene has the mutations).

--adapter_file 				To trim your short reads you need to specify what adaptors where used when sequencing. Arguments are  `TruSeq2-SE.fa`,
							 `TruSeq2-PE.fa`, `TruSeq3-PE.fa`. The default is `TruSeq3-PE.fa`.  

--genes_interest 			Path to a folder that contains all genes of interest. The correct formatting is ONE gene per file in FASTA format. This
							 folder can have any name, just make sure that it doesn't contain spaces in the name. Do not store additional files in 
							 this folder.

--assembly 					If you already have an assembly you can set this parameter and the pipeline will skip all the steps of assembly.

--run_classification  		The taxonomical classifications is based on a database of ~90G. This step will take several hours to complete downloading 
							the database. The default argument is `TRUE`. When `--run_classification FALSE` is not set, ProkGenomics will download 
							the database for taxonomical classification and stored in its base-directory. All runs after that will search in this 
							location for the database avoiding the lengthy download again. The download of the database depends on the available 
							disk space, this step will be skipped if there is not enough disk space.

--db_gtdb_path 				If you already download the gtdbtk database, indicate the path with this option. This parameter is only necessary if 
							you have the gtdbtk database in a different location from where ProkGenomics search by default (ProkGenomics/DB/db_gtdb/).

--keep_intermediate_files 	Keep intermediate files, including mapping files.

--cleanup 					The pipeline cleans-up by defaults. If you are debugging you can set the cleanup to FALSE and keep your work folders. 
							This will allow you to use `-resume` and troubleshoot errors.
"""
    exit 0
}


/*************************************************************************************************************************************************************
 * Channels
 read parameters in the command line 
 ************************************************************************************************************************************************************/

	Channel.value( params.cleanup )
		.set{ cleanup}

	Channel.value( params.run_classification )
		.set{ run_classification}

	Channel.value( params.keep_intermediate_files )
		.set{ keep_intermediate_files}

	//if(keep_intermediate_files){ params.keep_intermediate_files='copy'}else{params.keep_intermediate_files=null}

	Channel.value( params.threads )
		.set{ threads_ch }

	Channel.value( params.sample_name)
		.set{ sample_name}
	
	Channel.value( params.db_gtdb_path)
		.set{db_gtdb_path }

	Channel.value( params.percentage )
		.set{ percentage}

	shortreads = "${params.sample_path}/${params.sample_name}*{1,2}*"

	phage          ='phage'
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
if genes of interest a provided
*/
genes_fastas="${params.genes_interest}/*fasta"

if (params.genes_interest){ 
	Channel.fromPath( genes_fastas, checkIfExists: true )
                      .map { it -> [it.baseName, it] }
					  .set { ch_genes_interest }
	}


if (params.genes_interest){ 
	Channel.fromPath( genes_fastas, checkIfExists: true )
					  .collect()
                      .map { it -> [it.baseName, it] }
					  .set { ch_genes_interest_loop }
	}



//ch_genes_interest.view()

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
touch myDefaultInputFile_stats
touch myDefaultInputFile_qualimap
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

Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_stats")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_stats}


Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_qualimap")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_qualimap}


Channel
 	.fromPath("${params.default_empty_path}/myDefaultInputFile_assembly2gene")
	.map { file -> tuple (file.baseName,file) }
	.set{myDefaultInputFile_assembly2gene}
	

/**************************************************************************************************************************************************************
 * Print parameters
 **************************************************************************************************************************************************************/

//https://www.asciiart.eu/text-to-ascii-art
//ANSI SHADOW
log.info """
·········································································································
:  ██████╗ ██████╗  ██████╗ ██╗  ██╗ ██████╗ ███████╗███╗   ██╗ ██████╗ ███╗   ███╗██╗ ██████╗███████╗  :
:  ██╔══██╗██╔══██╗██╔═══██╗██║ ██╔╝██╔════╝ ██╔════╝████╗  ██║██╔═══██╗████╗ ████║██║██╔════╝██╔════╝  :
:  ██████╔╝██████╔╝██║   ██║█████╔╝ ██║  ███╗█████╗  ██╔██╗ ██║██║   ██║██╔████╔██║██║██║     ███████╗  :
:  ██╔═══╝ ██╔══██╗██║   ██║██╔═██╗ ██║   ██║██╔══╝  ██║╚██╗██║██║   ██║██║╚██╔╝██║██║██║     ╚════██║  :
:  ██║     ██║  ██║╚██████╔╝██║  ██╗╚██████╔╝███████╗██║ ╚████║╚██████╔╝██║ ╚═╝ ██║██║╚██████╗███████║  :
:  ╚═╝     ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═══╝ ╚═════╝ ╚═╝     ╚═╝╚═╝ ╚═════╝╚══════╝  :
·········································································································
ProkGenomics: Prokaryotic Genomics Pipeline
Author: Laura Perlaza-Jimenez (PhD) ~ Monash Genomics and Bioinformatics Platform
Rhys Grinter Laboratory
=========================================================================================================
	Type of assembly   : ${params.assembly_type}
	Short reads        : ${shortreads} 
	Long reads         : ${params.longreads}
	Assembly           : ${params.assembly}
	Reference          : ${params.reference}
	Output directory   : ${params.outdir}
	Number of threats  : ${params.threads}
	version            : ${params.version}       
==========================================================================================================
"""


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
include { db_checkv_download } from params.modules
include { checkv } from params.modules
include { plasclass } from params.modules

// Taxonomy classification
include { db_gtdb_download } from params.modules
include { gtdb } from params.modules

// split assembly
include { split_assembly } from params.modules

// Annotation
include { prokka } from params.modules
include { db_pharokka_download } from params.modules
include { pharokka } from params.modules

// Comparative Genomics
include { snippy } from params.modules
include { reference_format } from params.modules
include { minimap2 } from params.modules
include { get_coverage } from params.modules
include { qualimap } from params.modules

// Charaterization of genome
include { prodigal } from params.modules
include { makeblastdb } from params.modules
include { blast } from params.modules
include { extract_seq } from params.modules
include { assembly2gene } from params.modules
include { assembly2gene_table } from params.modules


//Create a final report of outputs
include { multiqc } from params.modules
include { report } from params.modules

//Clean up
include { cleanup} from params.modules

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

workflow db_checkv_download_workflow{
	main:
		db_checkv_download()

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


workflow db_pharokka_download_workflow{
	main:
		db_pharokka_download()

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
		reference_format(ch_in_reference)
		snippy(ch_in_reference,reference_format.out.reference_genome_ready,trimmed_reads)
		minimap2(reference_format.out.reference_genome_ready,trimmed_reads)
		get_coverage(minimap2.out.minimap2_path)
		qualimap(get_coverage.out.sorted_path)
	emit:
		snippy_path= snippy.out.snippy_path
		minimap2_path= minimap2.out.minimap2_path
		stats_path=get_coverage.out.stats_path
		qualimap_path=qualimap.out.qualimap_results
}


workflow db_gtdb_download_workflow{
	main:
		db_gtdb_download()
	emit:
		gtdbtk_db_path=db_gtdb_download.out.gtdbtk_db_path
		db_diskspace_command=db_gtdb_download.out.db_diskspace_command
}

workflow gtdb_workflow{
	take:
		contigs
		gtdbtk_db_path
		db_diskspace_command
	main:
		gtdb(contigs,gtdbtk_db_path,db_diskspace_command)
		//gtdb.out.gtdbtk_summary.view()
	emit:
		gtdb_path=gtdb.out.gtdbtk_summary
// PoisonPill check. Channel closed. Groovy solution

}

workflow makeblastdb_workflow{
	take:
		ch_genes_interest
	main:
		makeblastdb(ch_genes_interest)
	emit:
		blastDB=makeblastdb.out.blastDB

}


workflow prodigal_chr_workflow{
	take:
		sequence
		chromosome
	main:
		prodigal(sequence,chromosome)
	emit:
		prodigal_chr_path=prodigal.out.nucleotide_fasta

}

workflow prodigal_plasmid_workflow{
	take:
		sequence
		plasmid
	main:
		prodigal(sequence,plasmid)
	emit:
		prodigal_plasmid_path=prodigal.out.nucleotide_fasta

}

workflow prodigal_phage_workflow{
	take:
		sequence
		phage
	main:
		prodigal(sequence,phage)
	emit:
		prodigal_phage_path=prodigal.out.nucleotide_fasta

}


workflow blast_chr_workflow{
	take:
		sequence
		gene_id
		percentage
	main:
		sequence
			.combine(gene_id)
			.set{seq_gene}

		////seq_gene.view()
		blast(seq_gene,percentage)
	
	emit:
		blast_chr_path=blast.out.blast_output

}

workflow blast_plasmid_workflow{
	take:
		sequence
		gene_id
		percentage
	main:
		sequence
			.combine(gene_id)
			.set{seq_gene}

		////seq_gene.view()
		blast(seq_gene,percentage)
	emit:
		blast_plasmid_path=blast.out.blast_output

}

workflow blast_phage_workflow{
	take:
		sequence
		gene_id
		percentage
	main:
		sequence
			.combine(gene_id)
			.set{seq_gene}

		////seq_gene.view()
		blast(seq_gene,percentage)
	emit:
		blast_phage_path=blast.out.blast_output

}




workflow assembly2gene_chr_workflow{
	take:
		sequence
		blast_output
		gene_fasta
	main:
		sequence
			.combine(blast_output)
			.set{seq_gene}

		//seq_gene.view()
		extract_seq(seq_gene,gene_fasta)
		//extract_seq.out.extract_seqs_align.view()
		assembly2gene(extract_seq.out.extract_seqs_align)
		//assembly2gene.out.assembly2gene_nt_path.view()
	emit:
		assembly2gene_nt_path=assembly2gene.out.assembly2gene_nt_path
		assembly2gene_aa_path=assembly2gene.out.assembly2gene_aa_path
		assembly2gene_table_path=assembly2gene.out.assembly2gene_table_path


}

workflow assembly2gene_plasmid_workflow{
	take:
		sequence
		blast_output
		gene_fasta
	main:
		sequence
			.combine(blast_output)
			.set{seq_gene}

		////seq_gene.view()
		extract_seq(seq_gene,gene_fasta)
		//extract_seq.out.extract_seqs_align.view()
		assembly2gene(extract_seq.out.extract_seqs_align)
		//assembly2gene.out.assembly2gene_nt_path.view()
	emit:
		assembly2gene_nt_path=assembly2gene.out.assembly2gene_nt_path
		assembly2gene_aa_path=assembly2gene.out.assembly2gene_aa_path
		assembly2gene_table_path=assembly2gene.out.assembly2gene_table_path
}

workflow assembly2gene_phage_workflow{
	take:
		sequence
		blast_output
		gene_fasta
	main:
		sequence
			.combine(blast_output)
			.set{seq_gene}

		//seq_gene.view()
		extract_seq(seq_gene,gene_fasta)
		//extract_seq.out.extract_seqs_align.view()
		assembly2gene(extract_seq.out.extract_seqs_align)
		//assembly2gene.out.assembly2gene_nt_path.view()
	emit:
		assembly2gene_nt_path=assembly2gene.out.assembly2gene_nt_path
		assembly2gene_aa_path=assembly2gene.out.assembly2gene_aa_path
		assembly2gene_table_path=assembly2gene.out.assembly2gene_table_path

}

workflow assembly2gene_table_workflow{
	take:
	 	tables
	main:
		assembly2gene_table(tables)
	emit:
		assembly2gene_summary_path=assembly2gene_table.out.assembly2gene_summary_path

}


/// Is going thorugh just one gene!
/*
 * Defines the pipeline input parameters (with a default value for each one).
 * Each of the following parameters can be specified as command line options.
 
params.query = "$baseDir/data/sample.fa"
params.db = "$baseDir/blast-db/pdb/tiny"
params.out = "result.txt"
params.chunkSize = 100

db_name = file(params.db).name
db_dir = file(params.db).parent

workflow blast_workflow{
    /*
     * Create a channel emitting the given query fasta file(s).
     * Split the file into chunks containing as many sequences as defined by the parameter 'chunkSize'.
     * Finally, assign the resulting channel to the variable 'ch_fasta'
     
    Channel
        .fromPath(params.query)
        .splitFasta(by: params.chunkSize, file:true)
        .set { ch_fasta }

    /*
     * Execute a BLAST job for each chunk emitted by the 'ch_fasta' channel
     * and emit the resulting BLAST matches.
     
    ch_hits = blast(ch_fasta, db_dir)

    /*
     * Each time a file emitted by the 'blast' process, an extract job is executed,
     * producing a file containing the matching sequences.
     
    ch_sequences = extract(ch_hits, db_dir)

    /*
     * Collect all the sequences files into a single file
     * and print the resulting file contents when complete.
     
    ch_sequences
        .collectFile(name: params.out)
        .view { file -> "matching sequences:\n ${file.text}" }
}

*/

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
		minimap2_path
		stats_path
		gtdb_path
		assembly2gene_table
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
		minimap2_path,
		stats_path,
		gtdb_path,
		assembly2gene_table
		)
		
	emit:
		final_report=report.out.report_path
}

/*
workflow multiqc_workflow{
	take:
	 folder_name
	main:	
	 multiqc(folder_name)
	emit:
		multiqc_report=multiqc.out.multiqc_report
}

*/





/************************************************************************************************************************************************************** 
 * Main workflow
**************************************************************************************************************************************************************/



workflow{

	main:
	
	if (params.run_classification ){db_gtdb_download_workflow()}
	db_checkv_download_workflow()
	db_pharokka_download_workflow()

	if (params.assembly == null){
		if (params.assembly_type=='short'){
			shortreads_QC_workflow()
			shortreads_trim_workflow()
			shortreads_assembly_workflow(shortreads_trim_workflow.out.trimmed_reads)
				if (params.reference ){ comparative_genomics_workflow(ch_in_reference,shortreads_trim_workflow.out.trimmed_reads)
										snippy_output = comparative_genomics_workflow.out.snippy_path
										minimap2_output = comparative_genomics_workflow.out.minimap2_path
										stats_output = comparative_genomics_workflow.out.stats_path
										qualimap_output =comparative_genomics_workflow.out.qualimap_path
										}else{ 
										snippy_output = myDefaultInputFile_SNV_detection
										minimap2_output = myDefaultInputFile_mapping
										stats_output = myDefaultInputFile_stats
										qualimap_output = myDefaultInputFile_qualimap
										}

		
			assembly_qc_workflow(shortreads_assembly_workflow.out.scaffolds_path)
			extrachr_workflow(shortreads_assembly_workflow.out.scaffolds)
			split_assembly_workflow(shortreads_assembly_workflow.out.scaffolds,extrachr_workflow.out.plasclass_tsv,extrachr_workflow.out.checkv_summary)
			prokka_chr_workflow(split_assembly_workflow.out.chromosome_path,chromosome)
			prokka_plasmids_workflow(split_assembly_workflow.out.plasmid_path,plasmid)
			prokka_scaffolds_workflow(shortreads_assembly_workflow.out.scaffolds,denovoassembly)
			pharokka_workflow(split_assembly_workflow.out.phage_path)

			if (params.run_classification ){
				gtdb_workflow(split_assembly_workflow.out.chromosome_path,db_gtdb_download_workflow.out.gtdbtk_db_path,db_gtdb_download_workflow.out.db_diskspace_command)
				gtdb_output = gtdb_workflow.out.gtdb_path.ifEmpty{myDefaultInputFile_chr_classification }
				}else{ gtdb_output = myDefaultInputFile_chr_classification }


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
			stats_output = myDefaultInputFile_stats
			qualimap_output = myDefaultInputFile_qualimap

			extrachr_workflow(ch_in_assembly)
			split_assembly_workflow(ch_in_assembly,extrachr_workflow.out.plasclass_tsv,extrachr_workflow.out.checkv_summary)
			prokka_chr_workflow(split_assembly_workflow.out.chromosome_path,chromosome)
			prokka_scaffolds_workflow(ch_in_assembly,denovoassembly)
			if(plit_assembly_workflow.out.plasmid_path){prokka_plasmids_workflow(split_assembly_workflow.out.plasmid_path,plasmid)}
			if(split_assembly_workflow.out.phage_path){pharokka_workflow(split_assembly_workflow.out.phage_path)}

			if (params.run_classification ){
				gtdb_workflow(split_assembly_workflow.out.chromosome_path,db_gtdb_download_workflow.out.gtdbtk_db_path,db_gtdb_download_workflow.out.db_diskspace_command)
				gtdb_output = gtdb_workflow.out.gtdb_path.ifEmpty{myDefaultInputFile_chr_classification }
				}else{ gtdb_output = myDefaultInputFile_chr_classification }


			chromosome_path_output = split_assembly_workflow.out.chromosome_path.ifEmpty{ myDefaultInputFile_chr_extraction }
			plasmid_path_output = split_assembly_workflow.out.plasmid_path.ifEmpty{myDefaultInputFile_plasmid_extraction }
			phage_path_output = split_assembly_workflow.out.phage_path.ifEmpty{myDefaultInputFile_phage_extraction }
			prokka_scaffolds_path_output = prokka_scaffolds_workflow.out.prokka_path.ifEmpty{myDefaultInputFile_assembly_annotation }
			prokka_chr_path = prokka_chr_workflow.out.prokka_path.ifEmpty{myDefaultInputFile_chr_annotation }
			prokka_plasmids_path = prokka_plasmids_workflow.out.prokka_path.ifEmpty{myDefaultInputFile_plasmid_annotation }
			pharokka_path = pharokka_workflow.out.pharokka_path.ifEmpty{myDefaultInputFile_phage_annotation }
			plasclass_output = extrachr_workflow.out.plasclass_tsv.ifEmpty{myDefaultInputFile_plasmid_annotation }
			checkv_output = extrachr_workflow.out.checkv_summary.ifEmpty{myDefaultInputFile_phage_classification }
		}
 
		if(params.genes_interest){
			makeblastdb_workflow(ch_genes_interest)
			prodigal_chr_workflow(split_assembly_workflow.out.chromosome_path,chromosome)
			prodigal_plasmid_workflow(split_assembly_workflow.out.plasmid_path,plasmid)
			prodigal_phage_workflow(split_assembly_workflow.out.phage_path,phage)
			
			blast_chr_workflow(prodigal_chr_workflow.out.prodigal_chr_path,makeblastdb_workflow.out.blastDB,percentage)
			blast_plasmid_workflow(prodigal_plasmid_workflow.out.prodigal_plasmid_path,makeblastdb_workflow.out.blastDB,percentage)
			blast_phage_workflow(prodigal_phage_workflow.out.prodigal_phage_path,makeblastdb_workflow.out.blastDB,percentage)

			assembly2gene_chr_workflow(prodigal_chr_workflow.out.prodigal_chr_path,blast_chr_workflow.out.blast_chr_path,ch_genes_interest)
			assembly2gene_plasmid_workflow(prodigal_plasmid_workflow.out.prodigal_plasmid_path,blast_plasmid_workflow.out.blast_plasmid_path,ch_genes_interest)
			assembly2gene_phage_workflow(prodigal_phage_workflow.out.prodigal_phage_path,blast_phage_workflow.out.blast_phage_path,ch_genes_interest )
		
	 		table_chr=assembly2gene_chr_workflow.out.assembly2gene_table_path
			table_plasmids=assembly2gene_plasmid_workflow.out.assembly2gene_table_path
			table_phage=assembly2gene_phage_workflow.out.assembly2gene_table_path

	 		nt_chr=assembly2gene_chr_workflow.out.assembly2gene_nt_path
			nt_plasmids=assembly2gene_plasmid_workflow.out.assembly2gene_nt_path
			nt_phage=assembly2gene_phage_workflow.out.assembly2gene_nt_path

			aa_chr=assembly2gene_chr_workflow.out.assembly2gene_aa_path
			aa_plasmids=assembly2gene_plasmid_workflow.out.assembly2gene_aa_path
			aa_phage=assembly2gene_phage_workflow.out.assembly2gene_aa_path

			// make a list from the tuple, so the process is not input one by one but all at a time
			table_seqs=table_chr.concat(table_plasmids,table_phage).flatMap().collect()
			assembly2gene_table_workflow(table_seqs)
			assembly2gene_table=assembly2gene_table_workflow.out.assembly2gene_summary_path


		}else{assembly2gene_table=myDefaultInputFile_assembly2gene}


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
					 qualimap_output,
					 stats_output,
					 gtdb_output,
        			 assembly2gene_table

					)



/*	
	multiqc_workflow( sample_name,
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
					 stats_output,
					 gtdb_output,
					qualimap_output

					

        				//assembly2gene_table,
        				//assembly2gene_aligments,
        				//assembly2gene_peptides
					)
*/
	
cleanup(report_workflow.out.final_report)

}


workflow.onComplete {
	
	println ( workflow.success ? "\nDone! see the report in ${params.outdir} for more details \n" : "Oops .. something went wrong" )
}





//\n//Taxonomical classifications ${db_gtdb_download_workflow.out.db_diskspace_val} 