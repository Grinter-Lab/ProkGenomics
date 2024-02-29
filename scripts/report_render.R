
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
 novoassembly_path_QC= args[19]
 prokka_denovo_path= args[20]
 prokka_chr_path= args[21]
 prokka_plasmid_path= args[22]
 pharokka_path= args[23]
 plasmid_class= args[24]
 phage_class= args[25]
 snippy_path= args[26]

cmd <- paste(commandArgs(), collapse=" ")
cat("How R was invoked:\n");
cat(cmd, "\n")

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



file_name=paste0(prefix_sample,"_ProkGenomics_report.html")

rmarkdown::render( input  = report_template, 
                   params=list(
                                version= version, 
                                github= github,
                                fastqc_r1= fastqc_r1,
                                fastqc_r2= fastqc_r2,
                                fastqc_trim_r1= fastqc_trim_r1,
                                fastqc_trim_r2= fastqc_trim_r2,
                                novoassembly_path= novoassembly_path,
                                chromosome_path= chromosome_path,
                                plasmid_path= plasmid_path,
                                phage_path= phage_path,
                                prokka_denovo_path= prokka_denovo_path,
                                prokka_chr_path= prokka_chr_path,
                                prokka_plasmid_path= prokka_plasmid_path,
                                pharokka_path= pharokka_path,
                                plasmid_class=plasmid_class,
                                phage_class= phage_class,
                                novoassembly_path_QC= novoassembly_path_QC,
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



