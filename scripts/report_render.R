#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
report_template=args[1]
params_outdir= args[2]

parameters<- c( prefix_sample= args[3] ,
                version= args[4],
                github= args[5],
                fastqc_r1= args[6],
                fastqc_r2= args[7],
                fastqc_trim_r1= args[8],
                fastqc_trim_r2= args[9],
                novoassembly_path= args[10],
                chromosome_path= args[11],
                plasmid_path= args[12],
                phage_path= args[13],
                novoassembly_path_QC= args[14],
                prokka_denovo_path= args[15],
                prokka_chr_path= args[16],
                prokka_plasmid_path= args[17],
                pharokka_path= args[18],
                plasmid_class= args[19],
                phage_class= args[20],
                snippy_path= args[21],
                minimap2_path= args[22],
                gtdb_path= args[23])


parameters<-as.data.frame(t(gsub(".+myDefaultInputFile.+","NO_APPLY",parameters)))
file_name=paste0(parameters$prefix_sample,"_ProkGenomics_report.html")

cmd <- paste(commandArgs(), collapse=" ")
cat("How R was invoked:\n");
cat(cmd, "\n")
cat(args[1])

#print(paste0("if [ -s ", params_outdir,"/fastqc/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/fastqc/software_details.txt; fi"))

fastqc_version=try(system(paste0("if [ -s ", params_outdir,"/fastqc/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/fastqc/software_details.txt; fi"),intern = TRUE),silent=TRUE)
trimmomatic_version=try(system(paste0("if [ -s ",params_outdir,"/trimmomatic/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/trimmomatic/software_details.txt; fi"),intern = TRUE),silent=TRUE)
unicycler_version=try(system(paste0("if [ -s ",params_outdir,"/unicycler/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/unicycler/software_details.txt; fi"),intern = TRUE),silent=TRUE)
prokka_version=try(system(paste0("if [ -s ",params_outdir,"/prokka/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/prokka/software_details.txt; fi"),intern = TRUE),silent=TRUE)
pharokka_version=try(system(paste0("if [ -s ",params_outdir,"/pharokka/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/pharokka/software_details.txt; fi"),intern = TRUE),silent=TRUE)
plasclass_version=try(system(paste0("if [ -s ",params_outdir,"/plasclass/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/plasclass/software_details.txt; fi "),intern = TRUE),silent=TRUE)
checkv_version=try(system(paste0("if [ -s ",params_outdir,"/checkv/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/checkv/software_details.txt; fi"),intern = TRUE),silent=TRUE)
checkm_version=try(system(paste0("if [ -s ",params_outdir,"/checkm/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/checkm/software_details.txt; fi"),intern = TRUE),silent=TRUE)
snippy_version=try(system(paste0("if [ -s ",params_outdir,"/snippy/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/snippy/software_details.txt; fi"),intern = TRUE),silent=TRUE)
minimap2_version=try(system(paste0("if [ -s ",params_outdir,"/minimap2/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/minimap2/software_details.txt; fi"),intern = TRUE),silent=TRUE)
gtdb_version=try(system(paste0("if [ -s ",params_outdir,"/GTDB/software_details.txt ]; then sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/GTDB/software_details.txt; fi"),intern = TRUE),silent=TRUE)



rmarkdown::render( input  = report_template, 
                   params=list(
                                version= parameters$version, 
                                github= parameters$github,
                                fastqc_r1= parameters$fastqc_r1,
                                fastqc_r2= parameters$fastqc_r2,
                                fastqc_trim_r1= parameters$fastqc_trim_r1,
                                fastqc_trim_r2= parameters$fastqc_trim_r2,
                                novoassembly_path= parameters$novoassembly_path,
                                chromosome_path= parameters$chromosome_path,
                                plasmid_path= parameters$plasmid_path,
                                phage_path= parameters$phage_path,
                                prokka_denovo_path= parameters$prokka_denovo_path,
                                prokka_chr_path= parameters$prokka_chr_path,
                                prokka_plasmid_path= parameters$prokka_plasmid_path,
                                pharokka_path= parameters$pharokka_path,
                                plasmid_class=parameters$plasmid_class,
                                phage_class= parameters$phage_class,
                                novoassembly_path_QC= parameters$novoassembly_path_QC,
                                gtdb_path=parameters$gtdb_path,
                                minimap2_path=parameters$gtdb_path,
                                snippy_path=parameters$snippy_path,
                                fastqc_version= fastqc_version,
                                trimmomatic_version=trimmomatic_version,
                                unicycler_version=unicycler_version,
                                snippy_version=snippy_version,
                                prokka_version=prokka_version,
                                pharokka_version=pharokka_version,
                                plasclass_version=plasclass_version,
                                checkv_version=checkv_version,
                                checkm_version=checkm_version,
                                gtdb_version=gtdb_version,
                                minimap2_version=minimap2_version

                                ),
                        output_dir = Sys.getenv("PWD"),
                        output_file = file_name
                    )



