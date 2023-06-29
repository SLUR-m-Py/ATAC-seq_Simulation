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