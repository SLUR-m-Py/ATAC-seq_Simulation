## Load in needed mods
import numpy as np

## Ftn for difference of means 
def difference_of_means(x, y):
    """
    Calculate the difference in the means of two datasets x and y. 
    Returns a scalar equal to mean(y) - mean(x).
    """
    ## Return the dif in means
    return np.mean(y) - np.mean(x)

## Ftn for difference of means 
def difference_of_medians(x, y):
    """
    Calculate the difference in the means of two datasets x and y. 
    Returns a scalar equal to mean(y) - mean(x).
    """
    ## Return the dif in means
    return np.median(y) - np.median(x)

def difference_of_variance(x, y):
    """
    Calculates the difference in variance between x and y.
    """
    ## Return difference of variances
    return np.std(y)**2 - np.std(x)**2

## Ftn for permutations
def permute(x, y):
    """
    Given two datasets, return randomly shuffled versions of the concatonated inputs.
    """
    ## concatenate the data
    new = np.concatenate([x, y])
    ## shuffle the data
    np.random.shuffle(new)
    ## return the permuted data sets:
    return new[:len(x)], new[len(x):]

## Tests a null hypothesis 
def test_null(x, y, statistic, iters=1000):
    """
    Given two datasets (x and y), test a null hypothesis via permutation test with a given statistic.
    
    Params:
    x, y -- ndarrays, the data
    statistic -- a function of x and y
    iters -- number of times to bootstrap
    
    Ouput:
    a numpy array containing the bootstrapped statistic
    """
    ## Initilizse boots
    boots = []
    # Conduct the bootstrap
    while len(boots) < iters:
        boots.append(statistic(*permute(x, y)))
    ## Return the boots
    return np.array(boots)

## Ftn for non perametric boot
def bootstrap(x,ftn,n=1000,**kwargs):
    """
    Bootstraps a statsitic -- defined by input ftn --  on x. 
    """
    ## Initilze boots
    boots = []
    ## Iterate thru boot straps, randomly choose from x
    while len(boots) < n:
        boots.append( ftn( np.random.choice(x,size=len(x)), **kwargs ) )
    ## Return boots 
    return boots

## Ftn for calculating p-values
def calculate_pvalue(boots,obs):
    """
    Return the portion of the bootstraps larger than the observed.
    """
    ## Calculate the boots larger than the obs
    return sum(boots>=obs)/len(boots)


## Ftn for calculating ci 
def calculate_ci(boots,alpha=5):
    """
    Calculate the upper,lower percentiles and mean of input bootstraps.
    """
    ## calculate the lower and upper ci 
    lower,upper,mean = round(np.percentile(boots, alpha),4), round(np.percentile(boots,100-alpha),4),round(np.mean(boots),4)
    ## Print the mean and CI
    print(f'Mean = {mean}; CI = [{lower},{upper}]')
    ## Return the mean, lower and upper
    return mean,lower,upper