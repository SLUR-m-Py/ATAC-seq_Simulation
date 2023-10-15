#!/usr/bin/env python
#SBATCH --job-name=bamtofpkm            ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --partition=tb                  ## Size of partitions
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=13              ## Number of tasks to be launched
#SBATCH --mail-type=ALL                 ## Email for all job alerts
#SBATCH --mail-user=croth@lanl.gov      ## Email to this address
################################################
##     Fragment Counting Library & Script     ##
################################################
"""
Written by:

Cullen Roth, Ph.D.

Postdoctoral Research Associate
Genomics and Bioanalytics (B-GEN)
Los Alamos National Laboratory
Los Alamos, NM 87545
croth@lanl.gov
"""
## Set description of this library and scirpt
description = 'Calculates fragments along contigous bins within input bam files, standardizing them to fpkm, and saves them as a bed file.'

## Functions used throughout the script and for path manipulation and file checking
## Write a ftn for printing given a condition
def ifprint(message,bool):
    """Prints message given the boolean state."""
    ## Print the message
    print(message) if bool else None
    return bool 

## Write ftn for returning the parent path
def parentpath(inbam):
    """Returns the parent path of an input bam file."""
    ## Retrun the joined path
    return '/'.join(inbam.split('/')[:-1])

## Load in if-file ftn from os.path 
from os.path import isfile, basename

## Bring in ftns from os
from os import remove, environ, getcwd, makedirs

## Ftn for checking if we have a bam file
def isbam(inbam):
    """Checks if the input file is a .bam file"""
    ## Return the bool
    return ((inbam.split('.')[-1] == 'bam') and isfile(inbam))

## Ftn for removing a file
def throwout(filepath):
    """Given a path removes the file and prints to a log (if not None)."""
    ##  Take out the trash and if in verbose mode, leak the warn message print what we are doing
    return remove(filepath) if isfile(filepath) else None 

## Write a ftn for checking the computing environment
def condacheck(condaenv):
    """Checks the conda default environment matches that set in fridge."""
    ## Return a boolean on the check
    return environ['CONDA_DEFAULT_ENV'] == condaenv

## Ftn for submitting commands to bash and os, as well as multipfocessing
## Load in subprocess
import subprocess

## Ftn for calling the fastp, bwa and other commands to shell
def submitcall(command):
    """Submits a command to the shell via subprocess."""
    ## Submit the fastp command
    return subprocess.call(command,shell=True)

## Bring in multiprocessing
from multiprocessing import Pool

## Write ftn for parallizing commands
def topool(commands,npools):
    """Multiprocess commadns given the Pool module."""
    ## Make a pool to parallalize the commands, 
    with Pool(npools) as parall: ## map the calls of the command and return exit codes
        exit_codes = parall.map(submitcall, commands) ## Map the input command to the sbumitcall function
    ## Return the exit codes
    return exit_codes

## Functions for pysam and bam file manipulation
## Cuts on the bam extension
def rmdotbam(inbam):
    """Splits an input text on the string .bam"""
    ## Return the split on the sting
    return inbam.split('.bam')[0]

## Write a ftn to return the bam name
def justbam(inbam):
    """Returns the basename of an input bam file."""
    ## Returns the basename
    return basename(inbam)

## Write a ftn for indexing a bam file
def index(inbam,threads):
    """Submitts a samtools index command to index the input bam file (inbam)."""
    ## Format and submit a samtools command to index a bam file
    return f'samtools index -@ {threads} {inbam}'

## Write a ftn for checking if index is there
def hasix(inbam):
    """Checks if the index file .bai of an input bam file exists."""
    ## Return the boolean
    return isfile(inbam+'.bai') or isfile(inbam+'.csi')

## Load in pysam 
import pysam

## Loading in bams
def loadbam(bampath):
    """Loads in a bam file using pysam."""
    ## Run the alignment file command and return the bam object in binary format
    return pysam.AlignmentFile(bampath, "rb") 

## For getting chromlist
def getchrlist(sam_file,mito=['chrM']):
    """Given a bam file object returns a list of chromosome names, removing those stored in list MITO."""
    ## Return the filter chromosome list
    return [c for c in sam_file.references if c not in mito]

## Load in numpy 
import numpy as np 

## Getting chrom size from sam file
def chromlen(sam_file,chrom):
    """Returns a chromosome length using a bamfile object from pysam given a chromosome name (for e.g. chr1)."""
    ## Return the length given the chromosome name
    return np.array(sam_file.lengths)[(np.array(sam_file.references)==chrom)].min()

## Load in pandas
import pandas as pd 

## Ftn for counting fpkm along a chromosome
def chrombin(sam_file,chrom,window):
    """Generates windowed regions along a chromosome using a sam file object from pysam.\nThe size of regions is set by WINDOW."""
    ## Gather the sequence length and reset bedcounts
    seqlen, bed_counts = chromlen(sam_file,chrom), []
    
    ## Assert the sequence length is greater than zero
    assert (seqlen > 0), "ERROR: The chromosome has a length of zero!"
    
    ## Set a range given the window size and sequence length
    for j in range(1, seqlen, window):
        
        ## Set the stop given the window size
        stop = j+window-1 if j+window-1 < seqlen else seqlen
        
        ## Append the counts of a given region via the get reads ftn
        bed_counts.append((chrom, j, stop))   

    ## Format into a dataframe and return 
    return pd.DataFrame(bed_counts,columns=['Chrom','Left','Right'])

## Ftn for checking if a a file is a bed file
def isbed(inpath):
    """Returns a boolean if the input path ends in with bed extension and exists."""
    ## Returns the check
    return (inpath.split('.')[-1] == 'bed') and isfile(inpath)

## Ftn for loading in bed file
def loadbed(inpath,delim='\t',bednames=['Chrom','Left','Right']):
    """Loads an input bed file on input path."""
    ## Check the file is a bed file
    assert isbed(inpath), "ERROR: The input file -- %s -- does not exist or is not a .bed file!"%inpath
    ## Returns a tsv file
    return pd.read_csv(inpath,sep=delim,names=bednames)

## Ftn for removing bed extension
def basenobed(inbed):
    """Returns the basename of a bed file with no .bed file extension."""
    ## Return the the split 
    return basename(inbed).split('.bed')[0]

## Set defult vriables
myexcludes, windowsize, decimalsplace, threads = ['chrM'], 10000, 4, 12

## Set help variables
b_help = "Paths to input BAM files to count."
W_help = "Size in bp used to construct genomic windows, tiled per chromosome from left to right, for calculating fpkm (default: %s)."%windowsize
X_help = "Names of chromosomes to exclude from analysis (default: %s)."%', '.join(myexcludes)
O_help = "The name of the output path for saving the genome-count file generated here. This path will be made if it does not exist."
N_help = "The name of the output genomic count file. If none is provided chr1.%s.bed or genomic.%s.bed is used (as an example) for chromosome 1 with %s bp window size."%(windowsize,windowsize,windowsize)
T_help = "Number of threads to use in calls of samtools index (default: %s)."%threads
D_help = "Decimal place to round fpkm counts to (default: %s)."%decimalsplace
B_help = "Path to an input bed file defining regions to count alignments rather than uniform bins."
R_help = "Flag to standardize genomic counts to reads per kilobase per million."
F_help = "Flag to standardize genomic counts to fragments per kilobase per million."
C_help = "Flag to round fpkm counts to nearest whole integer."
V_help = "Flag to run in verbose mode. Default behavior is false."
S_help = "Flag to save out the calculated library sizes."

## If the library is called as an executable
if __name__ == "__main__":
    ## Load in argparser
    import argparse

    ## Set parser
    parser = argparse.ArgumentParser(description=description)

    ## Add required arguments 
    parser.add_argument("-b", "--bam-files",    dest="b", required=True,  type=str,  help=b_help, nargs='+')
        
    ## Add optional arguments
    parser.add_argument("-W", "--binsize",      dest="W", required=False, type=int,  help=W_help, default=windowsize,    metavar='n')
    parser.add_argument("-X", "--exclude",      dest="X", required=False, type=list, help=X_help, default=myexcludes,    metavar='chrM', nargs='+')
    parser.add_argument("-O", "--out",          dest="O", required=False, type=str,  help=O_help, default=None,          metavar=getcwd())
    parser.add_argument("-T", "--threads",      dest="T", required=False, type=int,  help=T_help, default=threads,       metavar='n')
    parser.add_argument("-D", "--n-decimals",   dest="D", required=False, type=int,  help=D_help, default=decimalsplace, metavar='n')
    parser.add_argument("-N", "--file-name",    dest="N", required=False, type=str,  help=N_help, default=None)
    parser.add_argument("-B", "--bed-file",     dest="B", required=False, type=str,  help=B_help, default=None)

    ## Add boolean variables
    parser.add_argument("-r", "--rpkm",         dest="R", help=R_help, action='store_true')
    parser.add_argument("-f", "--fpkm",         dest="F", help=F_help, action='store_true')
    parser.add_argument("-c", "--ceil",         dest="C", help=C_help, action='store_true')
    parser.add_argument("-V", "--verbose",      dest="V", help=V_help, action='store_true')
    parser.add_argument("-S", "--save-libsize", dest="S", help=S_help, action='store_true')

    ## Parse the arguments
    args = parser.parse_args()

    ## Set required vars
    bampaths = args.b

    ## Set variables from arg parser, such as the input paths etc
    window_bp, exclusions, outpath, filename, threads, ndecimals, bedpath = args.W, args.X, args.O, args.N, args.T, args.D, args.B

    ## Set boolean flags from arg parser
    verbose, torpkm, tofpkm, toceil, savesize = args.V, args.R, args.F, args.C, args.S

    ## Check all the input bam files are bam files
    for b in bampaths: ## Check if the input bam file is a bam file
        assert isbam(b), "ERROR: The input bam file -- %s -- is not a bamfile on this path: %s!"%(justbam(b),parentpath(b))

        ## Index the bam file 
        (index(b,threads-1),ifprint('INFO: Indexing input bam file: %s.'%justbam(b)),verbose) if not hasix(b) else None

    ## Reset outpath
    outpath = outpath if outpath else getcwd()

    ## Make the outpath
    makedirs(outpath,exist_ok=True)

    ## Load in the path exsits utility
    from os.path import exists as pathexists

    ## Make sure the outpath exists
    assert pathexists(outpath), "ERROR: The output path could not be found!"

    ## Check that the last chracter of the outpath is a forward slach "/"
    outpath = outpath if outpath[-1] == '/' else outpath + '/'

    ## Set the output name
    outname = outpath + filename if filename else outpath 

    ## Set the output path 
    tempoutname = outname + '.%s.%s.counts.bed' if filename else outname + '%s.%s.counts.bed'

    ## Load in a bamfile
    abamfile = loadbam(bampaths[0])

    ## Load in the list of chromosomes, initilize savepaths
    autosomes, savepaths, outpaths, bedcoms = getchrlist(abamfile, mito=exclusions), [], [], []

    ## Set the positional and bam columns names
    bamcolumns, positional = [justbam(b) for b in bampaths], ['Chrom','Left','Right']

    ## Format the column names
    thecolnames = positional + bamcolumns

    ## Bring in the bed file of defined regions 
    regions_df = loadbed(bedpath) if bedpath else None 

    ## Reset the window bp if a bed path was provided
    window_bp = basenobed(bedpath) if bedpath else window_bp

    ## Iterate thru the list of chromosome sin autosomes to make genomic bins
    for chrom in autosomes:
        ## Print the chrom if in verbose mode
        ifprint(chrom,verbose)

        ## Make genomic bins if given regions was passed
        chrombed = regions_df[(regions_df.Chrom==chrom)] if bedpath else chrombin(abamfile,chrom,window_bp)

        ## Set the savepaths 
        chromsavepath, chromoutpath = rmdotbam(bampaths[0])+'.%s.%s.bed'%(chrom,window_bp), tempoutname%(chrom,window_bp)

        ## Write to file
        chrombed.to_csv(chromsavepath,sep='\t',index=False,header=False)

        ## Format the call to bedtools multicov        
        call_to_bedtools = 'bedtools multicov -p -bams %s -bed %s > %s'%(' '.join(bampaths),chromsavepath,chromoutpath)

        ## Print the command to the output
        ifprint('INFO: %s'%call_to_bedtools, verbose)

        ## Append the chromosome bed savepath, the counts for chromosome outpaht, and the call to bed tools 
        savepaths.append(chromsavepath)
        outpaths.append(chromoutpath)
        bedcoms.append(call_to_bedtools)

    ## Gather the number of multipools we will need
    npools = len(bedcoms)

    ## Check our work
    assert npools > 0, "ERROR: The number of pools in multi processing was zero!"

    ## Submit the bedtools commands to python multiprocess
    exit_codes = topool(bedcoms,npools=threads)

    ## Check that all the exit codes are zero
    assert sum(exit_codes) == 0, "ERROR: Some of the calls (%s) to bedtools did not exectue properly!"%sum(exit_codes)

    ## Load in and merge the bed files
    count_df = pd.concat([pd.read_csv(filepath,sep='\t',names = thecolnames) for filepath in outpaths]).reset_index(drop=True)

    ## Make the sample columns integers
    for b in bamcolumns:
        count_df[b] = count_df[b].apply(int)

    ## Print the head
    ifprint('INFO:\n%s'%count_df.head(),verbose)

    ## Print the shape of the datframe
    ifprint('INFO: The dataframe has %s rows and %s columns.'%count_df.shape,verbose)

    ## if we need to standardize and we are not using an input bed path 
    if (torpkm or tofpkm or toceil) and (not bedpath):

        ## Calc the number of rows and columns, set the kb and million var, and gather the libsizes
        (nrows,ncols), kb, mil, libsizes = count_df.shape, 1000, 1000000, np.array([loadbam(b).count() for b in bampaths])

        ## Transform the libsizes'
        standards = (kb / window_bp) * (mil/libsizes)

        ## Make a dictionary, make sure the columns listed here match those in the other dataframe
        libsizedict = dict(zip(bamcolumns,standards))

        ## standardize by fpkm 
        if tofpkm or toceil: 
            ## Format the option            
            option = 'nearest whole fragment' if toceil else 'nearest %s decimal'%ndecimals

            ## Print what we are doing fpkm or ceil standardization
            ifprint('INFO: Standardizing counts to fpkm and rounding to ' + option, verbose)

            ## Iterate thru the columns
            for b in bamcolumns: # Iterate thru the sample names devide by two and round up to whole fragment
                count_df[b] = np.ceil(count_df[b].values/2)

            ## Normalize the counts via our libsizedict
            stand_df = count_df[bamcolumns].mul(libsizedict)

            ## iterate thru the columns
            for b in bamcolumns: # Iterate thru the sample names and round up 4 decimals
                stand_df[b] = np.ceil(stand_df[b].values) if toceil else np.round(stand_df[b].values,ndecimals) 

        elif torpkm: ## Print we are on rpkm 
            ifprint('INFO: Standardizing counts to rpkm.',verbose)

            ## Normalize the counts via our libsizedict
            stand_df = count_df[bamcolumns].mul(libsizedict)

        else: ## otherwise we are in trouble
            print("ERROR: This should not have happened but it did!")

        ## Merge the standard dataframe and the count
        count_df = pd.concat([count_df[positional],stand_df],axis=1)

        ## Check our work
        assert (count_df.shape[0] == nrows) and (count_df.shape[1] == ncols), "ERROR: There was an unexpected error in merging our standardized dataframe and positional columns from the raw counts dataframe."

        ## Print the transformed dataframe
        ifprint('INFO:\n%s'%count_df.head(),verbose)

        if savesize: ## Save out the libzies
            ## Set the saveout name for the library sizes
            libsavepath = outname+'.library.sizes.bed' if filename else outname+'library.sizes.bed'

            ## Print what we are doing
            ifprint('INFO: Saving library sizes of input bam files to file: %s'%justbam(libsavepath),verbose)

            ## Save out the lib sizes
            pd.DataFrame(list(zip(bampaths,libsizes)), columns=['Bamfile','Libsize']).to_csv(libsavepath,sep='\t',index=False)

    else: ## Don't standardize, return raw counts!
        ifprint('WARNING: Read counts are raw and not rpkm or fpkm standardized.',verbose)
        
    ## Format the save path
    saveoutpath = tempoutname%('genomic', window_bp)+'.gz'

    ## Print we are saving out the dataframe
    ifprint('INFO: Saving out genomic counts to: %s'%saveoutpath, verbose)

    ## Saveout the new count dataframe
    count_df.to_csv(saveoutpath,index=False,sep='\t', compression='gzip')

    ## Remove temporary files
    [(throwout(f),ifprint('WARNING: Removing file: %s'%justbam(f),verbose)) for f in outpaths+savepaths]

    ## Print we are done
    ifprint("INFO: Finished :-D",verbose)
## End of file 
