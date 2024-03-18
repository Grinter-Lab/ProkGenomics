#	my $rcmd=" Rscript $results_folder/scripts/curatingAlignments.r  
#$output_ref_chr 
#$gene_tag 
#$gene_genome_chr_out_alignment_nt_fasta 
#$gene_genome_chr_out_alignment_AA_fasta 
#$out_alignment_description --slave";

########################################################################################################################################################################
## load packages
########################################################################################################################################################################

rm(list=ls()); 
Sys.setenv('R_MAX_VSIZE'=32000000000)
suppressMessages(suppressWarnings(library(msa)))
suppressMessages(suppressWarnings(library(reshape2)))
suppressMessages(suppressWarnings(library(seqinr)))
suppressMessages(suppressWarnings(library( Biostrings)))
suppressMessages(suppressWarnings(library(stringr)))
options(warn=-1)
options(stringsAsFactors = FALSE)


########################################################################################################################################################################
## read arguments
########################################################################################################################################################################

args = commandArgs(trailingOnly=TRUE)
file=args[1]
genename=args[2]
out1=args[3]
out2=args[4]
out3=args[5]
setwd(dirname(file))
#print(getwd())
				
########################################################################################################################################################################
##Functions
########################################################################################################################################################################

	
findPotentialStartsAndStops <- function(DNA_string)
  { positions =c()
   	types =c()
     if(length(positions)==0){
     codons            <- c('atg', 'taa', 'tag', 'tga','ATG', 'TAA', 'TAG', 'TGA')
     for (i in  1:length(codons))
     {
        codon <- codons[i]
        occurrences <- matchPattern(codon, DNA_string)
      codonpositions <- occurrences\@ranges\@start  
        numoccurrences <- length(codonpositions) 	
   	positions   <- append(positions,codonpositions, after=length(positions))
     types       <- append(types,rep(codon, numoccurrences), after=length(types))
 }
    indices <- order(positions)
     positions <- positions[indices]
     types <- types[indices]
     mylist <- list(positions,types)
 return(mylist)
 
     }
   }     

 
  findORFsinSeq <- function(sequence)
  {
     mylist <- findPotentialStartsAndStops(sequence)
     positions <- mylist[[1]]
     types <- mylist[[2]]
     orfstarts <- numeric()
     orfstops <- numeric()
     orflengths <- numeric()
     numpositions <- length(positions)
 
     if (numpositions >= 2)
     {
        for (i in 1:(numpositions-1))
        {
           posi <- positions[i]
           typei <- types[i]
           found <- 0
           while (found == 0)
           {
              for (j in (i+1):numpositions)
              {
                 posj  <- positions[j]
                 typej <- types[j]
                 posdiff <- posj - posi
                 posdiffmod3 <- posdiff %% 3
           
                 orflength <- posj - posi + 3
                 if ((typei == 'atg' || typei == 'ATG') && (typej == 'taa' || typej == 'tag' || typej == 'tga'||typej == 'TAA' || typej == 'TAG' || typej == 'TGA') && posdiffmod3 == 0)
                 {
                    numorfs <- length(orfstops)
                    usedstop <- -1
                    if (numorfs > 0)
                    {
                      for (k in 1:numorfs)
                      {
                          orfstopk <- orfstops[k]
                          if (orfstopk == (posj + 2)) { usedstop <- 1 }
                      }
                    }
                    if (usedstop == -1)
                    {
                       orfstarts <- append(orfstarts, posi, after=length(orfstarts))
                       orfstops <- append(orfstops, posj+2, after=length(orfstops)) 
                       orflengths <- append(orflengths, orflength, after=length(orflengths))
                    }
                    found <- 1
                    break
                 }
                 if (j == numpositions) { found <- 1 }
              }
           }
        }
     }
     indices <- order(orfstarts)
     orfstarts <- orfstarts[indices]
     orfstops <- orfstops[indices]
     orflengths <- numeric()
     numorfs <- length(orfstarts)
     for (i in 1:numorfs)
     {
        orfstart <- orfstarts[i]
        orfstop <- orfstops[i]
        orflength <- orfstop - orfstart + 1
        orflengths <- append(orflengths,orflength,after=length(orflengths))
     }
     mylist <- data.frame(orfstarts=orfstarts,orfstops= orfstops,orflengths= orflengths)
     return(mylist)
  }
  
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
 		return(invisible(msa(x)))
            message('Successfully executed')
       		 },
        error = function(e){
            message('Caught an error!')
		return(invisible(msa(x)))

        },
        warning = function(w){
           message('Caught an warning!')
            #print(w)
        },
        finally = {
            message('All done, quitting.')
      		  }
    		)    
		}


		mySequences <- readAAStringSet(file)
		namesunique=unique(mySequences\@ranges\@NAMES)
		mySequences <- mySequences[namesunique]
		
		myFirstAlignment=Alignment(mySequences)

		location=ifelse(grepl('plasmid',file),'plasmid','chromosome')


		nSeq=length(rownames(myFirstAlignment))
		#print(paste('number found',nSeq,sep=' '))
	
  if(nSeq==2){
  		#print('2 sequences')
      		
     X_Alignment=myFirstAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))
				
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
			
			myAlignment=myFirstAlignment
			Ncopies=1
    }
    
      k=0
 	 if(nSeq==3){
 		#print('3 sequences')
 	 sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
  
    		
    		xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))

		
     X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))
				
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
######	
		
			if(any(xx[,2]=='-') && any(xx[,3]=='-')){Ncopies=1}else{Ncopies=2}		
   
}




k=0
 if(nSeq==4){
 		#print('4 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}

		xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		

     X_Alignment=myAlignment


			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))
			
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
######	

			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-')){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		
			
}


k=0
 if(nSeq==5){
 		#print('5 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
    				xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[5]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		     X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))
			
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
######	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-')){Ncopies=1}else{Ncopies=4}	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-' )){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		

}

k=0
 if(nSeq==6){
 		#print('6 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
    				xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[5]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[6]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		     X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')

######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))

			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
######	
			
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-') && any(xx[,6]=='-')){Ncopies=1}else{Ncopies=5}
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-')){Ncopies=1}else{Ncopies=4}	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-' )){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		

}

k=0
 if(nSeq==7){
 		#print('7 sequences')
 
  sequences=length(mySequences\@ranges\@NAMES)
 	  xx <- as.data.frame(matrix(ncol=sequences, nrow=myFirstAlignment\@unmasked\@ranges\@width[1]))
 
	  names=matrix(ncol=sequences, nrow=1)
	 
    		for( i in mySequences\@ranges\@NAMES){ 
        	x=c(unlist(strsplit(toString(unmasked(myFirstAlignment)[[i]]),'')))
        	k=k+1
    		xx[,k]=x
        	names[,k]=i
    		}
    				xx\$consensus=apply(xx, 1, function(x) {
						y=unique(x[-1])	
    					y=y[y!='-']	
						if(length(y)==0){y='-'}
						if(length(y)>1){if(x[1]==x[2]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[3]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[4]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[5]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[6]){y=x[1]}}
						if(length(y)>1){if(x[1]==x[7]){y=x[1]}}
						if(length(y)>1){y='N'}
						return(y)})

  		
  		conSeq=paste(xx\$consensus,collapse='')
		conSeq=gsub('-','',conSeq)
  		consensusSeq=DNAStringSet(conSeq)
  		consensusSeq\@ranges\@NAMES=mySequences\@ranges\@NAMES[2]
  		ref=DNAStringSet(mySequences[[1]])
  		ref\@ranges\@NAMES=mySequences\@ranges\@NAMES[1]
  		
  		myAlignment=invisible(Alignment(c(ref,consensusSeq)))
		     X_Alignment=myAlignment

			query=rownames(X_Alignment)[2]
      		Seq=toString(unmasked(X_Alignment)[[query]] )
      		Seq_string=strsplit(toString(Seq),'')
      		stops=which( Seq_string[[1]]=='-')
  		

			######
			ref=rownames(X_Alignment)[1]
      		Seqref=toString(unmasked(X_Alignment)[[ref]] )
      		Seqref_string=strsplit(toString(Seqref),'')
      		stopsref=which( Seqref_string[[1]]=='-')
      		
######	
			Align=data.frame(Ref=c(unlist(Seqref_string)),Query=c(unlist(Seq_string)))
				Align\$position=c(1:dim(Align)[1])

			Align\$Ref=as.character(Align\$Ref)
			Align\$Query=as.character(Align\$Query)
			Align\$mutation=paste(Align\$Ref,Align\$position,Align\$Query,sep='')
			Align\$Change='identical'
			Align\$Change[ Align\$Query!=Align\$Ref ]='SNP'
			Align\$Change[ Align\$Query=='-' ]='GAP'
			Align\$Change[ Align\$Ref=='-' ]='INSERTION'
			gaps=length(which(Align\$Query=='-'))
			
			changesSNPs=Align[Align\$Change=='SNP',]
			changesGAPs=Align[Align\$Change=='GAP',]
			changesINSERTION=Align[Align\$Change=='INSERTION',]
			SNPs=paste(changesSNPs\$mutation,collapse=',')
			GAPs=paste(changesGAPs\$mutation,collapse=',')
			INSERTION=paste(changesINSERTION\$mutation,collapse=',')
######	
		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-') && any(xx[,6]=='-') && any(xx[,6]=='-')){Ncopies=1}else{Ncopies=6}
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-') && any(xx[,6]=='-')){Ncopies=1}else{Ncopies=5}
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-') && any(xx[,5]=='-')){Ncopies=1}else{Ncopies=4}	
			if(any(xx[,2]=='-') && any(xx[,3]=='-') && any(xx[,4]=='-' )){Ncopies=1}else{Ncopies=3}		
			if(any(xx[,2]=='-') && any(xx[,3]=='-') ){Ncopies=1}else{Ncopies=2}		

}


names=mySequences\@ranges\@NAMES

if(Ncopies==1){

myAlignmentAA=invisible(msa(c(AAStringSet(Seqref),AAStringSet(Seq))))
RefAA=unlist(strsplit(toString(  toString(unmasked(myAlignmentAA)[[1]] ) ),''))
SeqAA=unlist(strsplit(toString(  toString(unmasked(myAlignmentAA)[[2]] ) ),''))
PepAlign=data.frame(Ref=unlist(RefAA),Query=unlist(SeqAA))
 
PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0

PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0

PepAlignSimilarity=sum(PepAlign\$match/length(PepAlign\$match))*100	

query=strsplit(names[2],'###')[[1]][1]
start=strsplit(names[2],'###')[[1]][2]
end=strsplit(names[2],'###')[[1]][3]
strand=strsplit(names[2],'###')[[1]][4]
GenomeFile=strsplit(names[2],':')[[1]][2]


startingcodon= min(which(PepAlign\$match==1))

if (stopsref != stops[1]){classification='possibly_truncated'}else{classification='present'}
if (startingcodon != 1){classification='late_startcodon'}else{classification='possibly_truncated'}


   		metadatasummary=data.frame(
      		gene=genename,
      			Ref=paste('\"',ref,'\"',sep=''),
 					Ref_length_nt=length(SeqAA)*3,
 					Ref_length_AA=length(RefAA),
 					Ref_stopCodon=stopsref,
 					GenomeFile=GenomeFile,
 					Query=paste('\"',query,'\"',sep=''),
 					Start=start,
 					End=end,
 					Strand= strand,
 					Query_length_nt=NA,
 					Query_length_AA=length(SeqAA),
 					nt_gaps=gaps,
 					stopCodon= stops[1],
 					SNPs=SNPs,
 					GAPs=GAPs,
 					INSERTION=INSERTION,
 					location=location,
 					classification=classification,
 					PepAlignSimilarity=PepAlignSimilarity,
 					DNAAlignSimilarity=NA,
 						Ncopies=Ncopies)

SeqQueryGenomeName=paste(GenomeFile,query,sep='::')
      		seq2Fasta(Seq,SeqQueryGenomeName,out1)
       		seqAA2Fasta(Seq,SeqQueryGenomeName, out2)
  write.table(metadatasummary,file=out3,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
}

count=1

if(Ncopies>1){

names=names[-1]

for (i in names){

count=count+1

query=strsplit(i,'###')[[1]][1]
start=strsplit(i,'###')[[1]][2]
end=strsplit(i,'###')[[1]][3]
strand=strsplit(i,'###')[[1]][4]
GenomeFile=strsplit(names[1],':')[[1]][2]

SeqName=xx[count]
ff=gsub('\"','',SeqName)
ff=gsub('\\n','',ff)
ff=gsub(', ','',ff)
ff=gsub('c\\\\(','',ff)
ff=gsub('\\\\)','',ff)

SeqName=ff
stops=which( Seq_trl_string[[1]]=='-')


myAlignmentAA=invisible(msa(c(AAStringSet(Seqref),AAStringSet(SeqName))))
RefAA=strsplit(toString(  toString(unmasked(myAlignmentAA)[[1]] ) ),'')
SeqAA=strsplit(toString(  toString(unmasked(myAlignmentAA)[[2]] ) ),'')
PepAlign=data.frame(Ref=unlist(RefAA),Query=unlist(SeqAA))
 
 
PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0
 
PepAlign\$match[PepAlign\$Ref== PepAlign\$Query]=1
PepAlign\$match[PepAlign\$Ref!= PepAlign\$Query]=0


PepAlignSimilarity=sum(PepAlign\$match/length(PepAlign\$match))*100	

startingcodon= min(which(PepAlign\$match==1))

if (stopsref != stops[1]){classification='possibly_truncated'}else{classification='present'}
if (startingcodon != 1){classification='late_startcodon'}else{classification='present'}

  		metadatasummary=data.frame(
      		gene=genename,
      			Ref=paste('\"',ref,'\"',sep=''),
 					Ref_length_nt=length(SeqAA)*3,
 					Ref_length_AA=length(RefAA),
 					Ref_stopCodon=stopsref,
 					GenomeFile=GenomeFile,
 					Query=paste('\"',query,'\"',sep=''),
 					Start=start,
 					End=end,
 					Strand= strand,
 					Query_length_nt=NA,
 					Query_length_AA=length(SeqAA),
 					nt_gaps=gaps,
 					stopCodon= stops[1],
 					SNPs=SNPs,
 					GAPs=GAPs,
 					INSERTION=INSERTION,
 					location=location,
 					classification=classification,
 					PepAlignSimilarity=PepAlignSimilarity,
 					DNAAlignSimilarity=NA,
 						Ncopies=Ncopies)

SeqQueryGenomeName=paste(GenomeFile,query,sep='::')
      		seq2Fasta(SeqName,SeqQueryGenomeName,out1)
       		seqAA2Fasta(SeqName,SeqQueryGenomeName, out2)
  write.table(metadatasummary,file=out3,append=TRUE,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
}


#print (genename)
#print (out3)
}
";