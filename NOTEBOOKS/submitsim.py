#!/usr/bin/env python

#SBATCH --job-name=submitsim            ## Name of job
#SBATCH --output=%x.%j.out              ## Name stdout
#SBATCH --error=%x.%j.err               ## Name stderr
#SBATCH --partition=tb                  ## Size of partitions
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
#SBATCH --mail-type=ALL                 ## Email for all job alerts
#SBATCH --mail-user=croth@lanl.gov      ## Email to this address

################################################
####    Submitting ATAC-seq Simulation      ####
################################################

## Bring in needed mods
import numpy as np, subprocess, pandas as pd

## Bring in mkdir and exsits ftns from os
from os import makedirs
from os.path import exists as pathexists

## Load in sleepy time
from time import sleep as sleepytime

## Ftn for calling the commands to shell
def submitcall(command):
    """Submits a command to the shell via subprocess."""
    ## Submit the input command
    return subprocess.call(command,shell=True)

## Ftn for calling subprocess
def subcheck(inputstr):
    """Checks output from subprocess submission of input str. Returns list of utf-8 coding."""
    ## call sub process
    return subprocess.check_output(inputstr, shell=True).decode("utf-8").split('\n')[:-1]

## Ftn for returning items in list greater than one
def biggerthanone(inlist):
    """Filters on items in list returing those with length larger than one."""
    ## Returns items in list
    return [k for k in inlist if len(k)>0]

## Make a slrum datarame
def slurmframe(user):
    """Generates a parsable dataframe from slurm queue for a given user."""
    ## Get the list of slurm user
    k = [biggerthanone(l.split(' ')) for l in subcheck('squeue -u %s'%user)]   
    ## Return the dataframe
    return pd.DataFrame(k[1:],columns=k[0])

## filter slurm df by job name and parition
def slurmfilter(user,jobname,partition):
    """Filters slurm fields given input user name on job id and partition."""
    ## Gather the slurm df
    slurm = slurmframe(user)
    ## Return the shape of the dataframe
    return slurm[(slurm.NAME==jobname) & (slurm.PARTITION==partition)].shape[0]

## Write ftn for making a directory
def dirmaker(dirpath):
    """Makes a directory given an input path."""
    ## Return the os command
    return makedirs(dirpath,exist_ok=True)

## Write ftn for checking file is a real deal!
def fileexists(filepath):
    """Checks the existance and size of an input file (foodstuffs)."""
    ## Via os, check if the path to "food" exists and meets our size threshold
    return pathexists(filepath)

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

## Write ftn for making count name
def countname(inbam,vp,rp,w,g,s):
    """Returns a file name for a genome count file given an input bam file, varying portion of peaks and reads sub-sampled."""
    ## Return the formated name
    return '%s_%s_%s_%s_%s_%s.tsv.gz'%(basebam(inbam).split('.bam')[0],vp,rp,w,g,s)

## Set node list
nodelist = ['gpu','fast','tb','mpi']

## jobs ber node grouop
jobpernode = 10

## Set the paritions size
jobsizes = jobpernode*np.ones(len(nodelist))

## Write ftn for submitting batchs to slurm
def topool(ftns,com='atac.sim',user='croth',wait=30,waitevery=200):
    """Submitss commands in batchs or locally."""    
    ## Initilize variables, the initizle comstatus, i, and partiton
    comstat, myiter, part, batchsize = 0, 0, [nodelist[-1]], sum(jobsizes)
    ## Iterativle add commands
    while (myiter < len(ftns)): ## So long as we have yet to submit all commands
        ## If the comstat is less than the batch size
        if (comstat < batchsize) and (len(part) > 0):
            ## Append the submission 
            submitcall('sbatch --partition %s '%part[0] + ftns[myiter]) 
            ## Add to to i
            myiter += 1
            if (((myiter+1) % waitevery) == 0) and ((myiter+1) < len(ftns)):
                while (comstat > 0):
                    ## Gather the current stats
                    comstat = sum([slurmfilter(user,com,n) for n in nodelist])
                ## sleep the script for a longer wait time 
                sleepytime(wait)
            else: ## Sleep for a bit
                sleepytime(wait/10)
        else: ## do nothing 
            pass 
        # Gather part counts per node 
        partcounts = [slurmfilter(user,com,n) for n in nodelist]
        ## Set the partitons
        part = [n for j,n in enumerate(nodelist) if partcounts[j] < jobsizes[j]]
        ## Redefine com status 
        comstat = sum(partcounts)
    ## Return subs 
    return myiter

## If the library is called as an executable
if __name__ == "__main__":

    ## ATAC trail 1 variables
    ## Set bam fle
    inputbam0 = '/panfs/biopan04/4DGENOMESEQ/ATAC_TRIAL1/2501_001/aligned/merged.2501_001.filtq30.bam'
    inputbam1 = '/panfs/biopan04/4DGENOMESEQ/ATAC_TRIAL1/2501_002/aligned/merged.2501_002.filtq30.bam'
    inputbam2 = '/panfs/biopan04/4DGENOMESEQ/ATAC_TRIAL1/2501_003/aligned/merged.2501_003.filtq30.bam'

    ## Set IDR bed
    inputbeda = '/panfs/biopan04/4DGENOMESEQ/ATAC_TRIAL1/ATAC_TRIAL1_chipr_peaks_all.bed'

    ## ----------------------------------------------- ##
    ## ATAC cryo trial 1 variables
    ## Set bam file
    inputbam3 = '/panfs/biopan04/4DGENOMESEQ/Cryo_ATAC/2501_007/aligned/merged.2501_007.filtq30.bam'
    inputbam4 = '/panfs/biopan04/4DGENOMESEQ/Cryo_ATAC/2501_008/aligned/merged.2501_008.filtq30.bam'

    ## Set IDR bed
    inputbedb = '/panfs/biopan04/4DGENOMESEQ/Cryo_ATAC/Cryo_TRIAL1_chipr_peaks_all.bed'

    ## ----------------------------------------------- ##
    ## ATAC cryo trial 2 variables
    ## Set bam file
    inputbam5 = '/panfs/biopan04/4DGENOMESEQ/SIMULATION/2501_018/aligned/merged.2501_018.filtq30.bam'
    inputbam6 = '/panfs/biopan04/4DGENOMESEQ/SIMULATION/2501_019/aligned/merged.2501_019.filtq30.bam'
    inputbam7 = '/panfs/biopan04/4DGENOMESEQ/SIMULATION/2501_020/aligned/merged.2501_020.filtq30.bam'

    ## Set IDR bed
    inputbedc = '/panfs/biopan04/4DGENOMESEQ/SIMULATION/Cryo_TRIAL2_chipr_peaks_all.bed'

    ## Group by bed and bam
    datagroups = [(inputbam0,inputbeda), (inputbam1,inputbeda), (inputbam2,inputbeda), (inputbam3,inputbedb), (inputbam4,inputbedb), (inputbam5,inputbedc), (inputbam6,inputbedc), (inputbam7,inputbedc)]

    ## Set the portion of varying peaks
    peak_portions = [0.01] + list(np.round(np.arange(0.05,1,0.05),2))

    ## Set the variables of the simulation, the total repliceates, the seed min
    read_portion, binsize, peaksize, replicates, seedn, calls, exists = 0.50, 10000, 0, 15, 100, [], 0

    ## Set the call to the ATAC-seq script 
    call = './simulation.py %s %s %s %s --verbose --seed %s --binsize %s --peaksize %s'

    ## Make error log dir
    dirmaker('./errors')

    ## Make the sims dir
    dirmaker('./sims')

    ## Iterate over the sample files 
    for (inputbam,inputbed) in datagroups:
        for r in range(replicates):
            ## Iterate over peak poritons 
            for p in peak_portions:
                
                ## format the file name
                outfilename = './sims/' + countname(inputbam,p,read_portion,peaksize,binsize,seedn)
                
                ## If the file was alread made 
                if fileexists(outfilename):
                    exists += 1
                else: ## Otherwise append thcall 
                    calls.append(call%(inputbam,inputbed,p,read_portion,seedn,binsize,peaksize))
                
            ## Update seed with plus one, this keeps autocorrelated peaks across a sample
            seedn += 1
            
    ## Check our work
    assert len(datagroups)*replicates*len(peak_portions) == len(calls) + exists, "ERROR: Unable to format needed calls!"
    
    ## Submit the calls to bash
    assert topool(calls) == len(calls), "ERROR: Not all of the calls of our simulation were submitted!"
## End of file 
