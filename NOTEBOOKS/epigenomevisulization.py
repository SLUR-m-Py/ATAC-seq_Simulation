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

## Import seaborn 
import seaborn as sns 

## Load in mystats lib
import mystatslib as sims 

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
    """Plots a bi-variate scatter plot of x vs y."""
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

def corrgrid(df,samples=None,figsize=(10,10),fontsize=12,xy=(0,10,30),
             name_dict=None,ytitlep=0.035,xlabel = None,alpha=0.2,color='k',kb=None):
    """Formats a grid of subplots to plot correlations."""
    ## Format the smaples
    toplot = samples if samples else getsamples(df)
    
    ## Calc the number to plot
    nplot = len(toplot)
    
    ## Calculat the kb
    kb = kb if kb else int(df.Right.min()/1000)
    
    ## Set the xlabel 
    myxlabel =  'Whole-Fragments per kb per Mil. (%s kb bins)'
    xlabel = xlabel if xlabel else myxlabel%kb
    
    ## Call figure, make the facecolor w
    fig,ax = plt.subplots(nplot,nplot,figsize=figsize,sharex=True,sharey=True)
    fig.set_facecolor('w')
    
    ## Iterate thru the samples 
    for i,m1 in enumerate(toplot):
        for j,m2 in enumerate(toplot):
            ## Plot if i greater than j 
            if (i > j):
                ## Set the axis
                plt.sca(ax[i,j])
                ## Gather the values 
                x,y = df[m1].values, df[m2].values
                ## Plot the data
                plt.plot(x,y,'o',color=color,alpha=alpha)
                ## IF xy triple was given
                if xy:
                    ## Set the x and y ticks 
                    plt.xticks(np.arange(np.min(xy),np.max(xy),np.median(xy)),fontsize=fontsize)
                    plt.yticks(np.arange(np.min(xy),np.max(xy),np.median(xy)),fontsize=fontsize)
                    ## Set the xlimits 
                    plt.xlim(np.min(xy),np.max(xy))
                    plt.ylim(np.min(xy),np.max(xy))
                else: ## Otherwise pass 
                    pass 
            ## If the sample is the sample, annotate it
            elif (i==j):
                ## Gather the sample name 
                s = name_dict[m1] if name_dict else m1
                ## Set the axis 
                plt.sca(ax[i,j])
                ## Add the title 
                plt.title(s,y=ytitlep,fontsize=fontsize)
                ## Turn off the axis 
                plt.axis('off')
            else:
                ## Turn off the axis 
                plt.sca(ax[i,j]);plt.axis('off')
            
    ## Annotate the x and y axis labels 
    fig.text(x=0.25,y=0.05,s=xlabel,fontsize=fontsize)
    fig.text(x=0.05,y=0.25,s=xlabel,fontsize=fontsize,rotation=90)
    ## Return the figure and axes 
    return fig,ax
    
def corrtable(df,samples=None,corr=None,cozero=True,ax=None,figsize=(7,5.5),vmin=0,vmax=1,cbar=True,cmap='coolwarm'):
    """
    Calcualtes pair-wise correlation metrics between samples in DF. 
    Plots a heatmap and returns a table of values. 
    """
    ## Format the smaples if none were passed 
    toplot = samples if samples else getsamples(df)

    ## Make into a dataframe 
    tmpdf = pd.DataFrame(np.ones((len(toplot),len(toplot))),columns=toplot,index=toplot)

    ## Iterate thru the samples 
    for i,m1 in enumerate(toplot):
        for j,m2 in enumerate(toplot):
            ## Plot if i greater than j 
            if (i > j):
                ## Gather the values 
                x,y = df[m1].values, df[m2].values
                ## Plot the data
                if corr == 'pearson':
                    c = sims.copearson(x,y,cozero=cozero)
                elif corr == 'nmi':
                    c = sims.conmi(x,y,cozero=cozero)
                elif corr == 'spearman':
                    c = sims.cospearman(x,y,cozero=cozero)
                elif corr == 'rsquared':
                    c = sims.cospearman(x,y,cozero=cozero)**2
                else:
                    c = sims.conmi(x,y,cozero=cozero)
                ## Add into the tmpdf in the upper and mirror diag
                tmpdf.loc[m1,m2] = c
                tmpdf.loc[m2,m1] = c
    ## Plot the results
    ## If ax was given do nothing 
    if ax:
        pass 
    else: ## otherwise plot a figure setting facecolro to 'w'
        fig,ax = plt.subplots(1,1,figsize=figsize)
        fig.set_facecolor('w')
    ## Set the axis 
    plt.sca(ax)
    ## Plot a headmap
    sns.heatmap(tmpdf,vmin=vmin,vmax=vmax,cbar=cbar,cmap=cmap,annot=True)

    ## Return the axis and table 
    return tmpdf, ax