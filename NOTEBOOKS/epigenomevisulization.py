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
###    Comet Plot Protocol & Visulization    ###
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
## Load in pandas and numpy
import pandas as pd, numpy as np

## Load in matploitlib for plotting
from matplotlib import pyplot as plt

## Load in glob
from glob import glob 

## Ftn for setting wc path
def sortglob(wc):
    ## Return the sorted glob of the wild card path 
    return sorted(glob(wc))

## Ftn for loading log
def loadlog(inpath):
    """Loads in and returns a log file."""
    ## Open with with
    with open(inpath,'r') as infile:
        data = infile.readlines()
        
    ## Return the data
    return data

## Deftn for parsing a log
def readstats(inpath):
    """Splits an input log file on path"""
    ## Return the split log
    return [k[:-1].split('\t') for k in loadlog(inpath)]

## Ftn for getting sample name 
def setsample(inpath):
    return inpath.split('/')[-1].split('.')[1]

## Ftn for making a dataframe
def dataframe(inpath):
    """Loads in mapping stats file and formats into a dataframe."""
    ## Make a temporary dataframe
    temp = pd.DataFrame(readstats(inpath))
    
    ## Set the sample name
    samplename = setsample(inpath)
        
    ## Set column names
    temp.columns = ['Mapping','Count_%s'%samplename,'Percent_%s'%samplename]
    
    ## Set the mapping as the index
    temp.index = temp.Mapping
    
    ## return the dataframe
    return temp

## Ftn for formating dataframe 
def formatdfs(inpaths):
    """Formats a dataframe values making percentages."""
    ## Load in dataframes, and initilzse list
    dfs,formatted  = [dataframe(p) for p in inpaths], []
    
    ## Iterate thru the dataframes
    for d in dfs:
        ## Set the last column
        s = d.columns[-1].split('Percent_')[-1]
    
        ## Format the column
        #d[s] = [ ('%s ( %s %s ) '%(r[d.columns[-2]],r[d.columns[-1]],'%') if j > 0 else r[d.columns[-2]]) for j,(i,r) in enumerate(d.iterrows())]
        d[s] = [r[d.columns[-2]] for j,(i,r) in enumerate(d.iterrows())]
        
        ## Append a copy of thta datframe
        formatted.append(d[[d.columns[-1]]].copy())

    ## return a concatinated list
    return pd.concat(formatted,axis=1).T

## Ftn for getting loc
def extractloci(inpath):
    return [int(k.split(' ')[1]) for k in loadlog(inpath) if 'Unique loci' in k][-1]

## Ftn for getting Frip score
def extractfrip(inpath):
    tmp = loadlog(inpath)
    frip_ix = [i+1 for i,l in enumerate(tmp) if 'INFO: The frip score' in l][0]
    return float(tmp[frip_ix][:-1].split('\t')[-1])

## Ftn for getting the number of unique loci per file
def getloci(inpaths):
    """Parses a run log to return the unique loci mapped via macs2."""
    ## Return the loci counts
    return pd.DataFrame([( setsample(l) , extractloci(l), extractfrip(l)) for l in inpaths], columns = ['Replicate','Loci','FrIP'])

## Ftn for getting just the upper tri
def undermask(df):
    ## Set a matrix on ones 
    mask = np.ones(df.shape,dtype='bool')
    ## Set the upper triangle index to false 
    mask[np.triu_indices(len(df))] = False
    ## Concat the mask 
    temp = np.concatenate(df.mask(~mask).values)
    ## Return the upper empty 
    return temp[~(np.isnan(temp))]

def loadbdg(inpath,chrom,left,right):
    
    temp = pd.read_csv(inpath,sep='\t',names=['Chrom','Start','End','FpKM'])
    temp = temp[(temp.Chrom==chrom)]
    temp['Meanpos'] = temp[['Start','End']].T.mean()
    temp = temp[(temp.Meanpos>=left) & (temp.Meanpos<=right)]
    return temp

## Ftn for checking if an object is not none type
def notnone(obj):
    """Checks if an object is None type."""
    ## Return the boolean
    return obj is not None

## Ftn for getting maximum counts
def getmax(x,y):
    """Returns the maximum count between a pair of input vectors x and y."""
    ## Return the max
    return np.max([np.max(x),np.max(y)])

## Ftn for getting samples
def getsamples(bed):
    """Returns all the column names of an input bed dataframe, thrid from the left."""
    ## Return all the columns, thrid from the left
    return bed.columns[3:].tolist()

## Ftn for loadding dataframe
def loaddf(inpath,sep='\t'):
    """Load in and return a tab seperated file."""
    ## Return the loaded csv
    return pd.read_csv(inpath,sep=sep)

## Ftn for making a paired_counts dataframe
def paired_counts(x,y):
    """Groups the input x and y values and counts them by unique pairs. Assumes x and y are descrete."""
    ## check the lenth of x and y
    assert len(x) == len(y), "ERROR: x and y do not have the same length."

    ## Make x and y pairs in to a dataframe with a zero column
    df = pd.DataFrame([x,y,np.zeros(len(x))],index=['X','Y','Counts']).T
    
    ## Count the uniqeu pairs of x and y
    counts = df.groupby(['X','Y']).count().reset_index()
    
    ## Log10 transforme those counts
    counts['Log10Counts'] = np.log10(counts.Counts.values)
    
    ## Return the counts
    return counts

## Ftn for plotting the data
def cometplot(x,y,plotmod=10,figsize=(6,5),ax=None,cbar_label='log$_{10}$ ( WFpkm Counts )',fontsize=12,tick_delta=5,xymod=1,cmap="jet"):

    ## Return the by counts dataframe and get the maximum of x and y
    bycounts, max_fpkm = paired_counts(x,y), getmax(x,y)

    ## Check if we need to make a figure and axes object
    ax = ax if notnone(ax) else plt.subplots(1,1,figsize=figsize)[1]

    ## Set the axis
    plt.sca(ax)
    
    ## Plot a one to one line given the maximum 
    plt.plot([0,max_fpkm],[0,max_fpkm],color='grey',linestyle='--',alpha=0.1)

    ## Set x and y -axis labels
    plt.xlabel('WFpkm',fontsize=12);plt.ylabel('WFpkm',fontsize=12)

    ## Via scatter plot ftn, plot the by counts dataframe
    sc = plt.scatter(bycounts.X.values,bycounts.Y.values,c=bycounts.Log10Counts.values,s=plotmod*(bycounts.Log10Counts.values+1),cmap=cmap)

    ## Set the x and y axis limits
    plt.xlim(-xymod,max_fpkm+xymod);plt.ylim(-xymod,max_fpkm+xymod)

    ## Adjust the x and y ticks
    plt.xticks(np.arange(0,max_fpkm+tick_delta,tick_delta),fontsize=fontsize-2);plt.yticks(np.arange(0,max_fpkm+tick_delta,tick_delta),fontsize=fontsize-2)
    
    ## Add a coloc bar
    cax = plt.colorbar(sc) if notnone(cbar_label) else None

    ## Set the color br label
    cax.set_label(cbar_label, fontsize=fontsize) if notnone(cax) else None

    ## Return the ax
    return ax

## WRite a ftn for plotting genes
def plotgene(gdf,y=0,lw=0.5,exonmod=5,ymod=0.05,fs=8):
    
    y = round(y,1)
    
    ## Gather the gene
    gene = gdf.Gene.min()
    
    ## Gather the strand
    strand = gdf[(gdf.Feature=='gene')].Strand.min()

    ## Color the gene by strand orientation
    color = 'tab:red' if strand == '-' else 'k'
    
    ## Plot the gene body
    plt.hlines(y,gdf.Start.min(),gdf.End.max(),linewidth=lw,color=color)
    
    ## Plot the exons
    exons = gdf[(gdf.Feature=='exon')]
    
    ## Iterate thru the exons
    [plt.hlines(y,j.Start,j.End,linewidth=lw*exonmod,color=color) for i,j in exons.iterrows()]
        
    ## Annotate the gene
    if strand == '-':
        plt.text(x=gdf.End.max(), y = y, s = gene, color = color, va='center',ha='left',fontsize=fs)
        
    else:
        plt.text(x=gdf.Start.min(), y = y, s = gene, color = color, va='center',ha='right',fontsize=fs)
        
    
    pass

## End of file 