#!/usr/bin/env python

#SBATCH --job-name=stat.res             ## Name of job
#SBATCH --output=./errors/%x.%j.out     ## Name stdout
#SBATCH --error=./errors/%x.%j.err      ## Name stderr
#SBATCH --time=10-20:00:00              ## Max time of submission (D,H,M,S)
#SBATCH --nodes=1                       ## Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ## Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ## Number of tasks to be launched
#SBATCH --mail-type=ALL                 ## Email for all job alerts
#SBATCH --mail-user=croth@lanl.gov      ## Email to this address

################################################
##     Co-zero Removal Statistic Library      ## 
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

## Bring in mods
import glob, pandas as pd, numpy as np, scipy.stats as ss

## Bring in norm. mutual info
from sklearn.metrics import normalized_mutual_info_score as nmi

## Ftn for finding co-zero index
def cozeroix(x,y):
    """Calculates the index where both arrays x and y are equal to zero (ie co-zero) and returns the boolean where they are not both zero."""
    ## Return the co-zero index
    return ~((x==0) & (y==0))

## Ftn for filtering out co-zeros
def filtcozero(x,y):
    """Returns x and y removing co-zeros."""
    ## Return the modified x and y
    return x[cozeroix(x,y)],y[cozeroix(x,y)]

## Write a ftn for ranking data
def rankdata(data):
    """Rank data using the scypi stats rank function and its default settings."""
    ## Returns the rank of input data using the ss.rankdata ftn.
    return ss.rankdata(np.array(data),axis=0)

## Tj for Kendalls w
def Tj(x):
    """Calculats the Tj(x) for use in Kendalls W."""
    ## Gather the counts of unique values in x
    u,l = np.unique(x,return_counts=True)
    ## Return the sum of the cubed ti in the counts if they are greater than 1
    return np.sum([(ti**3) - ti for ti in l[l>1]])

## ftn for Kdendalss Wt
def kendalls_wt(data):
    """Calculates kendalls W."""
    ## Rank the data
    data = rankdata(data)
    ## Return the rows and columns of input data
    n,m = data.shape
    ## Calculate R and Tj
    R,T = np.sum(data.sum(axis=1)**2), np.sum([Tj(t) for t in data.T])
    ## Calcualte little r and l
    r,l = (12*R) - ((3* n * m**2)*((n+1)**2)), ((n*m**2)*((n**2) - 1)) - (m * T)
    ## Return r over l
    return r/l

## A ftn cal to kendalls W, removing cozeros
def cokendallW(x,y,cozero=False):
    """Calls kendalls W, removing cozeros in x and y if True."""
    ##  return kendelts if co-zero is true, remove 
    return kendalls_wt(np.array([*filtcozero(x,y)]).T) if cozero else kendalls_wt(np.array([x,y]).T) 

## A ftn call to kendalls T, removing co zeros
def cokendallT(x,y,cozero=False):
    """Calculates Kendall T, removing cozeros in x and y if True."""
    ## Return kentall T if co-zero is true, remove zeros      
    return ss.kendalltau(*filtcozero(x,y))[0] if cozero else ss.kendalltau(np.array(x,dtype=float),np.array(y,dtype=float))[0]

## Ftn for returning norm. mutual info. 
def conmi(x,y,cozero=False):
    """Calcluates the normalized mutual information score between vectors x and y.\nIf cozero is True, this function removes cozeros in x and y."""
    ## Return the norm. mutual info.
    return  nmi(*filtcozero(x,y)) if cozero else nmi(x,y)

## A ftn call to spearman, removing cozeros
def cospearman(x,y,cozero=False):
    """Calculates spearmans P, removing cozeros in x and y if True."""
    ## if co-zero is true, remove
    return ss.spearmanr(*filtcozero(x,y))[0] if cozero else ss.spearmanr(x,y)[0]

## Ftn for returning pearson 
def copearson(x,y,cozero=False):
    """Calcluates the Pearson R correlation coeff. between vectors x and y.\nIf cozero is True, this function removes cozeros in x and y."""
    ## Calculating pearsons r, removing cozeros if needed
    return ss.pearsonr(*filtcozero(x,y))[0] if cozero else ss.pearsonr(x,y)[0]

## Define ftns for correlation
def myrank(x):
    """Ranks inpute data X using minimum method."""
    ## Returns my ranking method
    return len(x) + 1 - ss.rankdata(x, method='min')

## Ftn for calculating savage scores
def savage_score(ranks):
    """Calculate the savave scores of input RANKS."""
    ## Calculate the inverse of the total ranks
    l = 1/np.arange(1,len(ranks)+1,1) 
    ## Return their sum
    return np.array([np.sum(l[r-1:]) for r in ranks])   

## Write a ftn call to the top down ftn
def topdown(x,y):
    """Calcualtes the top down score of X and Y."""
    ## Return the pearson-r of the save score
    return ss.pearsonr(savage_score(myrank(x)),savage_score(myrank(y)))

## Define the call to the top down correlation and removing zeros
def cotopdown(x,y,cozero=False):
    """Calculates Top Down correlation, removing cozeros in x and y if True."""
    ## Calc top-down with zero remove if co-zero is true, else just the top down
    return topdown(*filtcozero(x,y))[0] if cozero else topdown(x,y)[0]

## Wrte ftn for loading in the dataframe
def loaddf(thepath,sep='\t'):
    """Loads in a tsv dataframe."""
    ## Load in and return the dataframe
    return pd.read_csv(thepath,sep=sep)

## Write ftn for loading in the ceil values of a given simulations
def raise_the_roof(thepath):
    """Loads in a genomic count dataframe and returns the rounded (ceiling) fpkm values from simulations."""
    ## Load in dataframe
    df = loaddf(thepath)
    ## Return the rounded values
    return df['Ceil_1'].values,df['Ceil_2'].values

## Write ftn for cunducting stats analysis
def calcstats(thepath, bools = [False,True]):
    """Conducts statistical tests given a the path to genomic counts simulation."""
    ## Gather x and y
    x,y = raise_the_roof(thepath)
    ## Calculate the top down correlation
    top_down, top_down_nz = [cotopdown(x,y,cozero=b) for b in bools]
    ## Calculate the person
    pearson, pearson_nz   = [copearson(x,y,cozero=b) for b in bools]
    ## The R2
    r2, r2_nz = pearson**2, pearson_nz**2 
    ## Calclualte the spearman
    spearman, spearman_nz = [cospearman(x,y,cozero=b) for b in bools]
    ## Calculate the kendall T
    kendallT, kendallT_nz = [cokendallT(x,y,cozero=b) for b in bools]
    ## Calculate kendalls W
    kendallW, kendallW_nz = [cokendallW(x,y,cozero=b) for b in bools]
    ## Calculate mutual information
    norm_mi, norm_mi_nz   = [conmi(x,y,cozero=b) for b in bools]
    ## Return the values of each statistic
    return  top_down, top_down_nz, pearson, pearson_nz, r2, r2_nz, spearman, spearman_nz, kendallT, kendallT_nz, kendallW, kendallW_nz, norm_mi, norm_mi_nz

## Set ftn to sort wild card file
def sortglob(inputwildcard):
    """Returns a list of sorted files from input wild card."""
    ## Return the sorted wild list
    return sorted(glob.glob(inputwildcard))

## Ftn for making a unique list
def uniquelist(x):
    """Returns a list of a set from input x."""
    ## Return a list of a set
    return list(set(x))

## if name is main 
if __name__ == "__main__":

    ## Load in argparse
    import argparse

    ## Set escriptions 
    description  = 'Caluclates pair-wise correlation statistics for simulated ATAC-seq replicates given an input genomicly-binned tsv file.'

    ## Set parser 
    parser = argparse.ArgumentParser(description=description)

    ## Add required arguments
    parser.add_argument('-i', '--input', dest='i', metavar='./path/to/genomic-bin.tsv.gz', required=True)

    ## Parse the arguments
    args = parser.parse_args()

    ## set output name
    inputpath = args.i

    ## Load in os path exists 
    from os.path import exists as pathexists

    ## Check the input path is real
    assert pathexists(inputpath), "ERROR: The input file -- %s -- cannot be found"%inputpath

    ## Check input path 
    assert (inputpath.split('.')[-2] == 'tsv') & (inputpath.split('.')[-1] == 'gz'), "ERROR: The input file is not a zipped tzv file!"

    ## format save out path
    saveoutpath = inputpath.split('.tsv')[0] + '.stats.csv'

    ## Bring in warnings
    import warnings

    ## Kendalls T throws a runtimewarning, ignore it
    warnings.filterwarnings("ignore", category=RuntimeWarning) 

    ## Set varables
    number_of_stats = 14

    ## Set statnames
    stat_names = ['TopDown','TopDownNz','Pearson','PearsonNz','Rsquared','RsquaredNz','Spearman','SpearmanNz','Kendallt','KendalltNz','Kendallw','KendallwNz','NormMutualInfo','NormMutualInfoNz']

    ## Check our work 
    assert len(stat_names) == number_of_stats, "ERROR: We are missing some stats!"

    ## Set experiment column names
    exp_columns = ['VaryingPeaks','ReadPortion','PeakSize','WindowSize','Seed']

    ## Gather the unique experiment name
    exp_name = inputpath.split('merged.')[-1].split('.')[0]

    ## Load in the simulation setings of each count
    tempdf = pd.DataFrame([inputpath.split('.tsv.')[0].split('_')[-5:]],columns=exp_columns)

    ## Iterate thru the dataframe and correct the column type
    for i,c in enumerate(tempdf.columns):
        ## Correct the column to float or integer
        tempdf[c] = tempdf[c].apply(float) if i in [0,1] else tempdf[c].apply(int)

    ## Add the experiment names 
    tempdf['ExpName'] = exp_name

    ## Iterate thru the simulations and return results as a dataframe
    the_stats = pd.DataFrame([calcstats(inputpath)], columns=stat_names)

    ## Merge with the temp dataframe
    results = pd.concat([tempdf,the_stats],axis=1)

    ## Save out results
    results.to_csv(saveoutpath,index=False)  
## End of file 