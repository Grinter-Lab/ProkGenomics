########################################################################################################################################################################
## load packages
########################################################################################################################################################################

rm(list=ls()); 
Sys.setenv('R_MAX_VSIZE'=32000000000)
suppressMessages(suppressWarnings(library(msa)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(stringr)))
options(warn=-1)

########################################################################################################################################################################
## read arguments
########################################################################################################################################################################

                                                                                                                                                                                             
args = commandArgs(trailingOnly=TRUE)
file=args[1]
fasta_file=args[2]

setwd(dirname(file))

########################################################################################################################################################################
##Functions
########################################################################################################################################################################

seq2Fasta <- function(seq,name,filename) 
				{
  						sink(filename,append=TRUE)
    					cat(paste0('>', name),file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)
    					the.sequence <- toString(seq)
    					the.sequence=gsub('N','-',the.sequence)
    					cat(the.sequence,file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)  
  						sink(NULL)
				}
				
						seqAA2Fasta <- function(seq,name,filename) 
				{
  						sink(filename,append=TRUE)
    					cat(paste0('>', name),file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)
    					the.sequence <- toString(seq)
    					cat(the.sequence,file=filename,append=TRUE)
    					cat('\\n',file=filename,append=TRUE)  
  						sink(NULL)
				}


Alignment <- function(x){
    	tryCatch(
        expr = {
 		return(invisible(msa(x,gapExtension =150,gap=1500,method='Muscle')))
            message('Successfully executed')
       		 },
        error = function(e){
          #  message('Caught an error!')
		return(invisible(  msa(x)))

        },
        warning = function(w){
            message('Caught an warning!')
            print(w)
        },
        finally = {
            message('All done, quitting.')
      		  }
    		)    
		}
########################################################################################################################################################################
## Main
########################################################################################################################################################################

mySequences <- readDNAStringSet(file)
nSeq=length(mySequences@ranges@NAMES)



if(nSeq>1){

		myFirstAlignment=Alignment(c(mySequences[1],mySequences[2]))
		myFirstAlignment_reverse=Alignment(c(reverseComplement(mySequences[1]),mySequences[2]))
		

 		FA=msaConsensusSequence(myFirstAlignment)
 		FAR=msaConsensusSequence(myFirstAlignment_reverse)

		Consensus=str_count(FA, pattern = '\\\\?')
 		Consensus_reverse=str_count(FAR, pattern ='\\\\?')
	
		if(Consensus<Consensus_reverse){myFirstAlignment=msa(mySequences)}else{myFirstAlignment=msa(c(reverseComplement(mySequences[1]),mySequences[-1]))}
 
  msaPrettyPrint(myFirstAlignment, output='tex', showConsensus = 'none', askForOverwrite=FALSE, verbose=FALSE,
                   alFile = fasta_file )

}

