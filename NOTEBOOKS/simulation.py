#!/usr/bin/env python

#SBATCH --job-name=atac.sim             ## Name of job
#SBATCH --output=./errors/%x.%j.out     ## Name stdout
#SBATCH --error=./errors/%x.%j.err      ## Name stderr
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=4               ## Number of tasks to be launched
"""
© 2023. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, 
and to permit others to do so.
"""
################################################
##########    ATAC-seq Simulation     ##########
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
## Load in needed mods
import pysam, pandas as pd, numpy as np, os, subprocess

## Bring in path exsits from os and make dirs
from os.path import exists as pathexists

## Set the various error messages used thru out the script
rep_err        = "ERROR: The number of peaks in each synthetic replicate is way off!"
var_err        = "ERROR: The number of varying peaks was larger than expected!"
chr_err        = "ERROR: Unable to generate a list of chromosomes!"
no_file_err    = "ERROR: We could not find the file: %s"
not_narrow_err = "ERROR: The bed file is not smaller than the input narrow peak file: %s, %s"
no_rows_err    = "ERROR: The bed dataframe has zero rows!"
make_rep_err   = "ERROR: The reads of replicate %s were not collected."
no_count_err   = "ERROR: The generated genome-wide fpkm file - %s - could not be found!"
fileis_err     = "ERROR: The output file - %s - already exists!\nDelete or move this file and try again."
sam_sub_err    = "ERROR: Unable to call samtools view command."
no_reads_err   = "ERROR: No reads were deteced from input bam file: %s"

## Set the various warning messages used in the script 
kept_bams = "WARNING: We retained the generated bam files:\n%s\n%s\n"
lib_warn  = "WARNING: The percent of library removed in replicate %s is: %s"

## Write a ftn for printing given a condition
def ifprint(message,bool):
    """Prints message given the boolean state."""
    ## Print the message
    print(message) if bool else None
    return bool 

## Write ftn for testing if an object is not equal to None
def notnone(obj):
    """A function to test if the input object is not equal to None"""
    ## Return the boolean if the obj is not equal to none
    return obj is not None

## Write ftn for checking file is a real deal!
def fileexists(filepath):
    """Checks the existance and size of an input file (foodstuffs)."""
    ## Via os, check if the path to "food" exists and meets our size threshold
    return pathexists(filepath)

## Write ftn for making a directory
def dirmaker(dirpath):
    """Makes a directory given an input path."""
    ## Return the os command
    return os.makedirs(dirpath,exist_ok=True)

## Ftn for checking if we have a bam file
def isabam(inbam):
    """Checks if the input file is a .bam file"""
    ## Return the bool
    return inbam.split('.')[-1] == 'bam'

## Get basename for bam 
def basebam(inbam):
    """Splits an input bam file by its path taking the basename"""
    ## Asser the file end is bam
    assert isabam(inbam), "ERROR: This is not a bam file!"
    ## Return the split
    return inbam.split('/')[-1]

## Write ftn for loading in bed
def loadbed(path,columns=['Chrom','Start','End'],sep='\t',header=None):
    """Loads in a bed file and returns the first three columns, the chromosome, start, and end."""
    ## Load in bed file
    bed = pd.read_csv(path,sep=sep,header=header).T[:3].T
    ## Set column names
    bed.columns = columns
    ## Return the bed
    return bed

## Write ftn for setting replicate names
def repname(inbam,rep,vp,rp,w,s):
    """Returns the name of a bam file given input bam file, the replicate number, the varying portion of peaks, and the portion of reads sampled from peaks."""
    ## Format the bam file name
    return './%s_%s_%s_%s_%s_%s'%(rep,vp,rp,w,s,basebam(inbam))

## Writ ftn for making count name
def countname(inbam,vp,rp,w,g,s):
    """Returns a file name for a genome count file given an input bam file, varying portion of peaks and reads sub-sampled."""
    ## Return the formated name
    return '%s_%s_%s_%s_%s_%s.tsv.gz'%(basebam(inbam).split('.bam')[0],vp,rp,w,g,s)

## Wrt ftn for slimming down narrow peak file
def slimdown(inbed,coi=['Chrom','Start','End']):
    """Returns a slimmer bed file from input, returning only columns of interst (COI) and droping duplicate rows."""
    ## Return only the columns of interest
    return inbed[coi].drop_duplicates()

## Wrt ftn for filltering bed file on chromosome list
def filt_chrom(inbed,chrlist,chrom='Chrom'):
    """Filteres input bed file by chrlist using pandas isin ftn."""
    ## Return rows with the column chrom in chrlist
    return inbed[(inbed[chrom].isin(chrlist))].copy()

## Calculating fpkm
def fpkm_norm(x,size,window=10000,mil=1000000,kb=1000):
    """Standardizes values in x to fpkm given a library size and window size (10 kb)."""
    ## Calc norm and return
    return x*(kb / window) * (mil/size)

## Loading in bams
def loadbam(bampath):
    """Loads in a bam file using pysam."""
    ## Run the alignment file command and return the bam object in binary format
    return pysam.AlignmentFile(bampath, "rb")

## Getting chrom size from sam file
def chromlen(sam_file,chrom):
    """Returns a chromosome length using a bamfile object from pysam given a chromosome name (for e.g. chr1)."""
    ## Return the length given the chromosome name
    return np.array(sam_file.lengths)[(np.array(sam_file.references)==chrom)].min()

## Get reads in sam file within a region of interest
def getreads(sam_file,c,l,r):
    """Returns the set of unique read names within a given region defined by a chromosome (c) and its left (l) and right (r) bounds."""
    ## Return the set of read names
    return set([read.query_name for read in sam_file.fetch(c,l,r)])

## Ftn for counting fpkm along a chromosome
def chromcount(sam_file,chrom,window=10000):
    """Count the number of reads across in windowed regions along a chromosome using a bamfile object from pysam.\nThe size of regions is set by WINDOW."""
    ## Gather the sequence length and reset bedcounts
    seqlen, bed_counts = chromlen(sam_file,chrom), []
    ## Assert the sequence length is greater than zero
    assert (seqlen > 0), "Error in chromosome length!"
    ## Set a range given the window size and sequence length
    for j in range(1, seqlen, window):
        ## Set the stop given the window size
        stop = j+window-1 if j+window-1 < seqlen else seqlen
        ## Append the counts of a given region via the get reads ftn
        bed_counts.append((chrom, j, stop, len(getreads(sam_file,chrom, j, stop))))   
    ## Format into a dataframe and return 
    return pd.DataFrame(bed_counts,columns=['Chrom','Left','Right','Counts'])

## Make a genome count file for a given bam file
def genomecount(sam_file,chrlist,window=10000,round_size=4):
    """Calculate a dataframe with fpkm counts across the input genome in CHRLIST."""
    ## Count the fragments per chromosome given an input bamfile object
    df = pd.concat([chromcount(sam_file,c,window=window) for c in chrlist])
    ## Get the lib size
    libsize = sam_file.count()
    ## Calculate the fpkm
    df['Norm'] = fpkm_norm(df.Counts, libsize)
    ## Round the norm counts up to the neares the whole integers
    df['Ceil'] = np.array(np.ceil(df.Norm), dtype=int)
    ## Round the norm to the nearst 4 decimals
    df['Norm'] = np.round(df['Norm'].values, round_size)
    ## Return the dataframe
    return df

## Sample reads within a list
def readsample(reads,p):
    """Randomly returns (without replacement) a set of read names from list given an input portion."""
    ## Return the randome set of reads
    return np.random.choice(list(reads),size=int(np.round(len(reads)*p)),replace=False)

## Sub sample reads in a given bed
def subsample(sam_file,bed,p):
    """Given loci listed witin the input bed file (with columns Chrom, Start, and End positions) returns a poriton (P) of reads within that region.\nA bamfile object from pysam is used for quickly parsing read names."""
    ## Sub samples a unique set or reads names in a genomic region (defined in input bed file)
    return set(np.concatenate([readsample(getreads(sam_file,r.Chrom,r.Start,r.End),p) for i,r in bed.iterrows()]))

## Randomly split a bed file
def ran_bed(inbed,p=0.50):
    """Randomly returns a portion of input bed dataframe given a poriton of its index."""
    ## Return the random index
    return inbed.loc[sorted(np.random.choice(inbed.index.values,size=int(np.round(len(inbed.index.values)*p)),replace=False)),:]

## Write a ftn for checking the computing environment
def openfridge(drawer):
    """Checks the conda default environment matches that set in fridge."""
    ## Return a boolean on the check
    return os.environ['CONDA_DEFAULT_ENV'] == drawer

## Write a ftn for expanding bed file
def expandbed(bed,window):
    """Given an input window size, centers, then expands the start and end postions of an input bed file."""
    ## Make a copy
    inbed = bed.copy()
    ## Take the middle and half way point
    inbed['Middle'], half = np.round(inbed[['Start','End']].mean(axis=1),0), int(window/2)
    ## Calculate new start and end positions
    inbed['Start'], inbed['End'] = inbed.Middle-half, inbed.Middle+half
    ## Return the moddified dataframe
    return inbed[['Chrom','Start','End']] if window > 0 else bed 

## Gather fragment names
def fragnames(inbam):
    """Returns the set of read names from input bam file."""
    ## Return the set, these are unique fragments
    return set([r.query_name for r in inbam])

## Ftn to write to file
def writetofile(inpath,intxt):
    """Opens a file to write lines to file."""
    ## With statment to open file
    with open(inpath,'w') as ofile:
        ofile.writelines(intxt)
    ## Return the path
    return inpath

## Format a samtools view command
def samtoolsview(i,o,opts='-Shb',threads=4):
    """Formats a samtools view command given input and output bam files, samtools view options and list of chromosomes."""
    ## Format and return the samtools command. See https://broadinstitute.github.io/picard/explain-flags.html for more details
    return f'samtools view -@ {threads} {opts} {i} > {o};samtools index -@ {threads} {o}'

## Ftn for calling the commands to shell
def submitcall(command):
    """Submits a command to the shell via subprocess."""
    ## Submit the input command
    return subprocess.call(command,shell=True)

## Ftn for adding '/' to end of path
def fixpath(inpath):
    """Appends a forward slash to an input path if missing from path."""
    ## fix and return the path
    return inpath if inpath[-1] == '/' else inpath + '/'

## Ftn for making a .txt file
def txtpath(inbam):
    """Formats an ouput txt file given an input bam file."""
    ## Format and return the new path and file name
    return inbam.split('.bam')[0]+'.txt'

## Ftn for formating columns
def formatcols(incolumns,rep):
    """Reformats columns given the list of columns and replicate number (1 or 2)."""
    ## Return the reformated clumns
    return incolumns[:3] + ['%s_%s'%(c.split('_')[0],rep) for c in  incolumns[3:]]

def cleanup(toremove):
    """Removes files listed in toremove"""
    ## Return the removed file status
    return [os.remove(k) for k in toremove]

## Set needed variables
description  = 'Run simulation of ATAC-seq replicates given input bam and narrow peak file.'
myscriptsdir = '/panfs/biopan04/4DGENOMESEQ/EPICPIPELINE'
libsize_tol  = 0.5
ranseed      = 117 
hicenv       = 'hicexplorerenv'
peaksize     = 0 
windowsize   = 10000
samthreads   = 4

## If the library is called as an executable
if __name__ == "__main__":

    ## Load in argparser
    import argparse

    ## Set parser
    parser = argparse.ArgumentParser(description=description)

    ## Add required arguments
    parser.add_argument('bam',                     metavar='./path/to/input/file.bam',    type=str,   help='A path to an input bam file used to construct synthetic replicates.')
    parser.add_argument('bed',                     metavar='./path/to/input/narrowPeak',  type=str,   help='A path to an input narrowPeak (bed) file from macs2 used to construct synthetic replicates.')
    parser.add_argument('varying_peaks_portion',   metavar='V',                           type=float, help='The portion (0 - 1) of peaks from input bed files to vary between replicates.')
    parser.add_argument('reads_sampled_portion',   metavar='P',                           type=float, help='The portion (0 - 1) of reads from peaks to be sampled out of selected peaks.')
    
    ## Add optional arguments
    parser.add_argument("--seed",                  metavar='n',                           type=int,   help="Sets the random seed (default: %s) to be used in numpy."%ranseed,                                                                                             default=ranseed)
    parser.add_argument("--out",                   metavar='./output/path',               type=str,   help="The name of the output path for saving the genome-count file generated here.",                                                                                default='./sims')
    parser.add_argument("--env",                   metavar='hicexploererenv',             type=str,   help="Name of python environment for executing this script (default: %s)."%hicenv,                                                                                  default=hicenv)
    parser.add_argument("--peaksize",              metavar='n',                           type=int,   help="Size in bp of window formed around peak center (default: %s).\nPassing zero leaves input bed file unchanged, using peak bounds as the window size."%peaksize, default=peaksize)
    parser.add_argument("--binsize",               metavar='n',                           type=int,   help="Size in bp used to construct genomic windows, tiled per chromosome from left to right, for calculating fpkm (default: %s)."%windowsize,                       default=windowsize)
    parser.add_argument("--threads",               metavar='n',                           type=int,   help="Number of threads set in samtools view command (default: %s)"%samthreads,                                                                                     default=samthreads)

    ## Add boolean variables
    parser.add_argument("--verbose",  help="Flag to run in verbose mode. Default behavior is false.",                 action='store_true')
    parser.add_argument("--keepbams", help="Flag to keep bam files of synthetic replicates generated from analysis.", action='store_true')

    ## Parse the arguments
    args = parser.parse_args()

    ## Set variables from arg parser, set the input paths, the portion of varying peaks, reads sampled out of varying peaks, the window size for the counting of fpkm
    thebampath, narrowpath,  varying_peaks_portion, reads_sampled_portion, window_bp = args.bam, args.bed, args.varying_peaks_portion, args.reads_sampled_portion, args.binsize 

    ## Set verbosity and if we are keeping the bams made from simulation
    verbose, keepbams, theenv = args.verbose, args.keepbams, args.env

    ## Check the environment is set
    openfridge(theenv)

    ## Print the arguments if in verbose mode
    ifprint(args,verbose)

    ## Reset the exclusion list of chromosomes, the peak size if none was passed, the random seed
    peakwindow, seed, pathout  = args.peaksize, args.seed, args.out

    ## Set the seed
    np.random.seed(seed)

    ## make the dir
    dirmaker(pathout)

    ## Set the output path  
    outpath = fixpath(pathout) + countname(thebampath,varying_peaks_portion,reads_sampled_portion,peakwindow,window_bp,seed)

    ## Check that the output file in gzip or tsv form dosen't exist
    assert ((not fileexists(outpath)) and (not fileexists(outpath.split('.gz')[0]))), fileis_err%outpath

    ## Print the arguments
    [ifprint(str(j),verbose) for j in ['\nArguments passed to this script:\n',thebampath, narrowpath, varying_peaks_portion, reads_sampled_portion, verbose, seed]]

    ## Check that the bam ane bed files are real
    for f in [thebampath,narrowpath]:
        assert fileexists(f), no_file_err%f

    ## Load in the bamfile and calculdate its size
    bamfile = loadbam(thebampath)

    ## Set the bamfile 
    readnames = fragnames(bamfile)

    ## Make sure there are reads
    assert len(readnames) > 0, no_reads_err%thebampath

    ## Print number of reads
    ifprint(len(readnames),verbose)

    ## Calculate the libsize 
    libsize = bamfile.count()

    ## print the lib size
    ifprint(libsize,verbose)

    ## Set the chromosome list taking only autosomes
    autosomes = ['chr%s'%(i+1) for i in range(22)]

    ## Check our work
    assert len(autosomes) >= 0, chr_err

    ## Load in the narrow peak file
    narrow = loadbed(narrowpath)

    ## Slim down the narrow peak file
    slimbed = slimdown(narrow)

    ## Expand the bed 
    expbed = expandbed(slimbed,peakwindow)

    ## Filter the bed file on chromosome
    bed = filt_chrom(expbed,autosomes)

    ## Check that the bed file is smaller than the narrow peak file, check rows
    assert bed.shape[0] <= narrow.shape[0], not_narrow_err%(bed.shape[0],narrow.shape[0])

    ## Check that the bed file is smaller than the narrow peak file, check columns
    assert bed.shape[1] <= narrow.shape[1], not_narrow_err%(bed.shape[1],narrow.shape[1])

    ## Remove the bed files we no longer need in memory
    del expbed, slimbed, narrow 

    ## Gather the portion of peaks ot vary
    varying_peaks = ran_bed(bed,p=varying_peaks_portion)

    ## Split the varying peaks randomly into two sets
    r1_vary = ran_bed(varying_peaks)

    ## Take the not listed peaks
    r2_vary = varying_peaks[~(varying_peaks.index.isin(r1_vary.index))]

    ## Check the shape
    assert np.abs(r1_vary.shape[0] - r2_vary.shape[0]) < 3, rep_err

    ## Check our work
    assert r1_vary.shape[0] <= np.ceil(varying_peaks.shape[0]/2), var_err

    ## Check that peaks were selected
    assert np.min([r1_vary.shape[0],r2_vary.shape[0]]) > 0, no_rows_err

    ## Gather the sets of reads to remove to make synthetic replicates 1 and 2 
    set_k1, set_k2 = subsample(loadbam(thebampath),r1_vary,reads_sampled_portion), subsample(loadbam(thebampath),r2_vary,reads_sampled_portion)
    
    ## Print the len of set
    [ifprint(s,verbose) for s in (len(set_k1),len(set_k2))]

    ## What are the percentages of reads we are taking of the total library size?
    p_lib1, p_lib2 = len(set_k1)/libsize, len(set_k2)/libsize

    ## Round the size of the libraries
    round_p_lib = np.round([p_lib1,p_lib2],2)

    ## Print the percent lib size removed for each replicate made
    [ifprint(lib_warn%(i+1,j),verbose) for i,j in enumerate(round_p_lib)]

    ## Check that the percentages of the taken libraries are within tolerance 
    assert np.diff(round_p_lib) < libsize_tol

    ## Make the sets for each new synthetic repliceate, removing the varying individual read sets
    set1_readnames, set2_readnames = readnames - set_k1, readnames - set_k2

    ## Check our work, that the sets are not empty 
    assert len(set1_readnames) > 0, make_rep_err%1
    assert len(set2_readnames) > 0, make_rep_err%2

    ## Print the lenths of the total read sets
    [ifprint(s,verbose) for s in (len(set1_readnames),len(set2_readnames))]

    ## Set the replicate paths
    r1_path, r2_path = [repname(thebampath,r,varying_peaks_portion,reads_sampled_portion,peakwindow,seed) for r in ['r1','r2']]

    ## Set the names of read paths
    r1_txt_path, r2_txt_path = txtpath(r1_path), txtpath(r2_path)

    ## Write the read sets to file for replicate 1 and 2
    read1_path, read2_path = writetofile(r1_txt_path, [r+'\n' for r in list(set1_readnames)]), writetofile(r2_txt_path, [r+'\n' for r in list(set2_readnames)])

    ## Format the samtools commands to make the replicates
    read1_com, read2_com = samtoolsview(thebampath, r1_path, opts='-Shb -N %s'%read1_path, threads=4), samtoolsview(thebampath, r2_path, opts='-Shb -N %s'%read2_path, threads=4)

    ## Submit the samtools view commands to the os 
    assert submitcall(read1_com) + submitcall(read2_com) == 0, sam_sub_err

    ## Print the lenths of the total read sets
    [ifprint(s,verbose) for s in (read1_com,read2_com)]

    ## Make the genome count dataframe for rep 1 and rep 2
    r1_bed, r2_bed = genomecount(loadbam(r1_path), autosomes, window=window_bp),  genomecount(loadbam(r2_path), autosomes, window=window_bp)

    ## Reformat and set columns of each replicate
    r1_bed.columns, r2_bed.columns = formatcols(r1_bed.columns.tolist(),1), formatcols(r2_bed.columns.tolist(),2)

    ## merge the files togther
    genomebed = r1_bed.merge(r2_bed)

    ## Save out the file
    genomebed.to_csv(outpath,sep='\t',compression='gzip',index=False)

    ## Check that the output file is real
    assert fileexists(outpath), no_count_err%outpath

    ## See if we need to clean up, if so remove the large bam files
    k = ifprint(kept_bams%(r1_path,r2_path),verbose) if keepbams else cleanup([r1_path,r2_path,r1_path+'.bai',r2_path+'.bai',r1_txt_path,r2_txt_path])
## End of file 
