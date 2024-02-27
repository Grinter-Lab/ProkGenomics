
args = commandArgs()
 report_template=args[6]
 params_outdir= args[7]
 prefix_sample= args[8] 
 version= args[9]
 github= args[10]
 fastqc_r1= args[11]
 fastqc_r2= args[12]
 fastqc_trim_r1= args[13]
 fastqc_trim_r2= args[14]
 novoassembly_path= args[15]
 chromosome_path= args[16]
 plasmid_path= args[17]
 phage_path= args[18]
 prokka_denovo_path= args[19]
 prokka_chr_path= args[20]
 prokka_plasmids_path= args[21]
 pharokka_path= args[22]
 plasmid_class= args[23]
 phage_class= args[24]
#cmd <- paste(commandArgs(), collapse=" ")
#cat("How R was invoked:\n");
#cat(cmd, "\n")


tryCatch(fastqc_version=system(paste0("sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/fastqc/software_details.txt")),error=function(e) e)
tryCatch(trimmomatic_version=system(paste0("sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/trimmomatic/software_details.txt")),error=function(e) e)
tryCatch(unicycler_version=system(paste0("sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/unicycler/software_details.txt")),error=function(e) e)
tryCatch(prokka_version=system(paste0("sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/prokka/software_details.txt")),error=function(e) e)
tryCatch(pharokka_version=system(paste0("sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/pharokka/software_details.txt")),error=function(e) e)
tryCatch(plasclass_version=system(paste0("sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/plasclass/software_details.txt")),error=function(e) e)
tryCatch(checkv_version=system(paste0("sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/checkv/software_details.txt")),error=function(e) e)
tryCatch(checkm_version=system(paste0("sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/checkm/software_details.txt")),error=function(e) e)
tryCatch(snippy_version=system(paste0("sed -n '/version/, /version/{ /version/!p }' ", params_outdir,"/snippy/software_details.txt")),error=function(e) e)



file_name=paste0(prefix_sample,"_ProkGenomics_report.html")

          rmarkdown::render( input  = report_template, 
                                        params=list(
                                                    version= quote(version), 
                                                    github= quote(github),
                                                    fastqc_r1= quote(fastqc_r1),
                                                    fastqc_r2= quote(fastqc_r2),
                                                    fastqc_trim_r1= quote(fastqc_trim_r1),
                                                    fastqc_trim_r2= quote(fastqc_trim_r2),
                                                    novoassembly_path= quote(novoassembly_path),
                                                    chromosome_path= quote(chromosome_path),
                                                    plasmid_path= quote(plasmid_path),
                                                    phage_path= quote(phage_path),
                                                    prokka_denovo_path= quote(prokka_denovo_path),
                                                    prokka_chr_path= quote(prokka_chr_path),
                                                    prokka_plasmids_path= quote(prokka_plasmids_path),
                                                    pharokka_path= quote(pharokka_path),
                                                    plasmid_class=quote(plasmid_class),
                                                    phage_class= quote(phage_class),
                                                    fastqc_version= fastqc_version,
                                                    trimmomatic_version=trimmomatic_version,
                                                    unicycler_version=unicycler_version,
                                                    snippy_version=snippy_version,
                                                    prokka_version=prokka_version,
                                                    pharokka_version=pharokka_version,
                                                    plasclass_version=plasclass_version,
                                                    checkv_version=checkv_version,
                                                    checkm_version=checkm_version
                                                    ),
                                        output_dir = Sys.getenv("PWD"),
                                        output_file = file_name
                                        )



