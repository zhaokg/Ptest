from . import Rbeast as cb
from .cvt_to_numpy import force_convert_to_numpy  ### Y = force_convert_to_numpy(Y)

def beast123( Y,                  # Y can be 1D (one time series), 2D (many ts), or 3D (stacked images)
              metadata = None,    # an object with many fields to describe Y 
              prior    = None,    # an object with many fields to specify the bayesian model
              mcmc     = None,    # an object to specify the MCMC simulation parameters 
              extra    = None,    # an object to specify the outputs 
              method   = 'bayes', # 'bayes','bic','aic','aicc','hic'
             ):
      """
   
   
################################################################################################
 Bayesian changepoint detection and time series decomposition for regular or irregular time series data
    
   The fitted model is:
	 Y= trend + error             if data has no periodic/seasonal variation (i.e., season='none')
     Y= trend + seasonal + error  if data has periodic/seasonal variation 
     Y= trend + outlier  + error  if data is trend-only (no seasonal variation) but with potential outliers
     Y= trend + seasonal + outlier + error if data has periodic/seasonal variation and also has outliers
   where trend is a piecewise linear or polynomial function with an unknown number of trend changepoints to 
   be inferred; seasonal is a piecewise periodic function with an unknown number of seasonal changepoints to 
   be inferred; and the outlier component refers to potential pikes or dips at isolated data points and is included
   only if metadata.hasOutlierCmpnt=True (in beast123) or hasOutlier=True ( in beast or beast_irreg)


*Note*:
--------------------------------------------------------------------------------------------------
    (1) beast123() is an all-inclusive function that duplicates the functionalities of beast() 
     and beast_irreg(). beast123() can handle a single, multiple, or 3D stacked time series, being
     either regular or irregular.  beast123() internally supports parallel computing if the input is mutliple time
     series or stacked time series images. It allows for customization through four object argumentsa: metadata,
     prior, mcmc, and extra.

    (2) Currently, the BEAST algorithm is implemented to handle regular time series purely for computational 
    advantages. The function 'beast_irreg' does accept irregular time series but internally it still 
    aggregates them into regular ones prior to applying the BEAST model. For the aggregation, both the "time"
    and "deltat" args are needed to specify indvidial times of data points and the regular time interval desired. 

*Quick Examples*:
--------------------------------------------------------------------------------------------------
import Rbeast as rb
    
nile,yr =  rb.load_example('nile')             # annual flow of the Nile river
rb.beast_irreg( nile, time=yr, season='none' ) # time should be explicitly provided for beast_irrg
    
ohio = rb.load_example('ohio')                # a satelite ndvi time series at an Ohio site
o    = rb.beast_irreg(ohio.ndvi, time=ohio.rdate, deltat=1/12, period=1.0)
rb.plot(o) 
    
    
*Input arguments*:
--------------------------------------------------------------------------------------------------
    Y:  an iregularly-spaced time series or 1D signal; it should be a numeric vector. For reggular 
    time series, use 'beast' or 'beast123' instead. For multiple time series  or stacked time series 
    images such as satellite data, use 'beast123'.
 
     ... :  the remaining arguments are many paired keywords and values to specifiy time information 
    or parameters for the beast algorithm. Check the R version of BEAST for detailed explanations
    (https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf). Below is a brief description.
 
*Possible Keywords*:
--------------------------------------------------------------------------------------------------   
 
     |> metadata <| a struct variable specifying metadata for the input Y 

     metadata.season           : a string specifier about the seasonal component. Possible values are
	                             (1) 'none' - trend-only data with no seasonlity or periodic variation,
                                 (2) 'harmonic' - harmonic model for the seasonal component,     
                                 (3) 'dummy' - dummy model for the seasonal component     
                                 (4) 'svd': (experimental feature) the same as 'harmonic' except that the 
                                  periodic/seasonal component is modeled as a linear combination of function 
                                  bases derived from a Single-value decomposition. The SVD-based basis functions 
                                  are more parsimonious than the harmonic sin/cos bases in parameterizing the 
                                  seasonal variations; therefore, more subtle changepoints are likely to be detected.
     metadata.period           : the period of the seasonal/periodic component
     metadata.isRegularOrdered : Deprecated. Now use 'metadata.startTime' and ' metadata.deltaTime' to 
                                 specify times for evenly-spaced data, and use ' metadata.time' to
                                 provide times of individual data points for irregular or regular data	 
     metadata.startTime        : the start time of regular input
     metadata.deltaTime        : the time interval between consecutive data
                                 points. For regular time seires, it
                                 determins the times of all the data points;
                                 for irregular time series, deltaTime
                                 specifies the interval at which the
                                 irregular input is aggregated/resampled
     metadata.time             : for irregular inputs (or regular or unordered data), specify 
                                 the individual times of data points. Bothstring and numeric vectors
                                 are allowed, check rdrr.io/cran/Rbeast/man/beast123.html
                                 for more details.
     metadata.whichDimIsTime   : which dim of the 2D or 3D input refer to time
     metadata.missingValue     : a value indicating bad/missing data
     metadata.maxMissingRate   : the max missingness rate beyond which the ts will be ignored
     metadata.deseasonalize    : if true, the input ts is first de-seasonalize by removing a global 
                                seasonal component, prior to applying BEAST
     metadata.detrend          : if true, the input ts is first de-trended by removing
                                 a global trend component, prior to applying BEAST
     metadata.hasOutlierCmpnt  : boolean; if true, the model with an outlier component.
                                 Y = trend + outlier + error if season='none' or
                                 Y = trend+season+outlier+error will be fitted.
								 
 
    |> prior <|  prior hyperparameter for the BEAST model
  
     prior.seasonMinOrder  : min harmonic order considered for the seasonal compnt
     prior.seasonMaxOrder  : max harmonic order considered for the seasonal compnt
     prior.seasonMinKnotNum : min number of seasonal changepoints allowed
     prior.seasonMaxKnotNum : max number of seasonal changepoints allowed
     prior.seasonMinSepDist : the min segment length of a seaosnal segment
                             (i.e, the min seperation distance between
                             neighorboring seasonal changepoints)
     prior.seasonLeftMargin:  an integer;  the number of leftmost data points excluded for 
	                         seasonal changepoint detection. That is,  no seasonal changepoints are allowed in 
	                         the starting window/segment of length seasonLeftMargin. seasonLeftMargin must
	                         be an unitless integer–the number of time intervals/data points so that the
	                         time window in the original unit is seasonLeftMargin*deltat. If missing,
	                         seasonLeftMargin defaults to the minimum segment length 'seasonMinSepDist'
     prior.seasonRightMargin:  an integer;  the number of rightmost data points excluded for 
	                         seasonal changepoint detection. That is,  no changepoints are allowed in 
	                         the ending window/segment of length seasonRightMargin. seasonRightMargin must
	                         be an unitless integer–the number of time intervals/data points so that the
	                         time window in the original unit is seasonRightMargin*deltat. If missing,
	                         seasonRightMargin defaults to the minimum segment length 'seasonMinSepDist'
     prior.trendMinOrder	: min polyonimal order considered for the trend compnt
     prior.trendMaxOrder	: max polyonimal order considered for the trend compnt
     prior.trendMinKnotNum  : min number of trend changepints allowed
     prior.trendMaxKnotNum  : max number of trend changepints allowed
     prior.trendMinSepDist  : the min segment length of a trend segment(i.e, the min seperation distance
                             between any two neighorboring trend changepoints)
     prior.trendLeftMargin:  an integer;  the number of leftmost data points excluded for 
	                         trend changepoint detection. That is,  no trend changepoints are allowed in 
	                         the starting window/segment of length trendLeftMargin. trendLeftMargin must
	                         be an unitless integer–the number of time intervals/data points so that the
	                         time window in the original unit is trendLeftMargin*deltat. If missing,
	                         trendLeftMargin defaults to the minimum segment length 'trendMinSepDist'
     prior.trendRightMargin:  an integer;  the number of rightmost data points excluded for 
	                         trend changepoint detection. That is,  no trend changepoints are allowed in 
	                         the ending window/segment of length trendRightMargin. trendRightMargin must
	                         be an unitless integer–the number of time intervals/data points so that the
	                         time window in the original unit is trendRightMargin*deltat. If missing,
	                         trendRightMargin defaults to the minimum segment length 'trendMinSepDist'							 
    prior.outlierMaxKnotNum : an integer; needed only if meta.hasOutlierCmpnt=True to specify the maximum 
                             number of outliers (i.e., outlier-type changepoints) allowed in the time series	 

     |> mcmc <|    parameters to set up the MCMC sampler
 
     mcmc.seed             : the seed for the random number generator            
     mcmc.samples          : the number of samples to be collected
     mcmc.thinningFactor   : chain thinning factor; take every
                             'thinningFactor'-th sample from the chain
     mcmc.burnin           : number of initial samples discarded
     mcmc.chainNumber      : number of chains
     mcmc.maxMoveStepSize  : a moving step
     mcmc.trendResamplingOrderProb  : probability to re-sample trend order
     mcmc.seasonResamplingOrderProb : probability to re-sample seasonal order
     mcmc.credIntervalAlphaLevel    : signifiance level for computing
                                      credible interval
 
     |> extra <|    parameters to control computaiton or outputs

     extra.dumpInputData      : if true, dump the input time series (o.data)
     extra.whichOutputDimIsTime: which dim of the output array refer to time
     extra.computeCredible     : compute credible intervals (unless needed
                                 for assessing the curve fitting, it is
                                 better to disable it for changepoint
                                 analysis because computing CI is
                                 expensive.)(o.trend.CI and o.season.CI)
     extra.fastCIComputation   : if true, a fast version of algorithm used to
                                compute CI. But it is still expensive.
                                Set computeCredible=false if CI is of litte
                                interest
     extra.computeSeasonOrder  : if true, estimate harmonic orders needed to
                                adequately fit the seasonal component
                                (o.season.order)
     extra.computeTrendOrder   : if true, estimate polynomial orders needed to
                                adequately fit the trend component
                                (o.trend.order)
     extra.computeSeasonChngpt  : if true, dump the scp detected
                                (o.season.cp, o.season.cpCI, and cpAbruptChange)
     extra.computeTrendChngpt   : if true, dump the tcp detected
                                (o.trend.cp, o.trend.cpCI, and cpAbruptChange)
     extra.computeSeasonAmp     : if true, estimate the seasonal ammplitue
                                (o.seasin.amp)
     extra.computeTrendSlope    : if true, estimate time-varying slope and
                                 the associated probabilities of slopes being positive or negative
     extra.tallyPosNegSeasonJump:
     extra.tallyPosNegTrendJump :
     extra.tallyIncDecTrendJump :
     extra.printProgressBar     : if true, show a progress bar
     extra.printOptions         : if true, print the BEAST parameters at the
                                 start
     extra.consoleWidth         : the console/terminal width
     extra.numThreadsPerCPU     : IMPORTANT for remote sensing/earth science
                                 applicaitons: specify the numbers of
                                 concurrent threads per cpu core
     extra.numParThreads        : specify the number of total concurrent
                                threads
     |> method <|  a string specifier of the method for formulating model posterior probability.
          Possible values are 
         (1) 'bayes': the full Bayesian formulation (this is the default)  
         (2)'bic': approximation of posterior probability using the Bayesian information criterion (bic)
         (3)'aic': approximation of posterior probability using the Akaike information criterion (aic)
         (4)'aicc': approximation of posterior probability using the corrected Akaike information criterion (aicc)
         (5)'hic': approximation of  posterior probability using the Hannan–Quinn information criterion  (hic)
         
		 
     *Note*:
    --------------------------------------------------------------------------------------------------
     beast, beast_irreg, and beast123 calls the same dll extension library written in C/C++
    (i.e., Rbeast.pyd); they are just wrappers to the core BEAST algorithm.
    In particular, BEAST currently handles only regular time series;
    irregular inputs will be first aggregated  before applying BEAST.

 
    *Result/Output*: The output is a struct variable; example of the fields include
    --------------------------------------------------------------------------------------------------
        marg_lik: marginal likilood; the larger, the better
        sig2    : variance  of error
        trend   : the trend component; a struct variable (say, T)
        season  : the season componet; a stuct variable  (say,S)
        The subfields of trend or season:
        .ncpPr        : the prob distribution for number of changepoints
        .ncp          : mean number of changepoints in trend or seasonality
        .ncp_meidan   : median number of changepoints
        .ncp_mode     : mode from ncpPr
        .ncp_pct90    : 90# percentile from ncpPr
        .cpOccPr      : changepoint occurrance probability over time
        .cp           : list of all possible changepoints (many are not sigficant)
        .cpPr         : occurrence probability of the changepoints in cp
        .cpAbruptChange: the sudden changes in trend or seasonlity at cp
        .cpCI         : confidence interval of the cps
        .Y            : the fitted trend or seasonality 
        .SD           : standard deviation of the fitted Y
        .CI           : Credible interval of the fittted Y
        .order   : the mean harmonic or polynomial orders estimated to fit the seasonal and trend     
        trend.slp     : slope of the trend 
        trend.slpSD   : standard dev of the estimated slope
        trend.slpSgnPosPr: time-varying probability of the slope being postive
        trend.slpSgnZeroPr: time-varying probability of the slope being 0
        season.amp     : amplitue of the estiamted seasonality overtime
        season.ampSD   : standard ev of the estiamated amplitude
 
More help:  
--------------------------------------------------------------------------------------------------
       This terse help doc sucks (I know); so far, the best details are still the
       R help doc, available at https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf.
       Python doesn't allow a '.' in variable names, so Python's equivalent to R's 
       beast(Y,start=1987,tcp.minmax=c(0,5)) is beast(Y, start=1987,  tcp_minmax=[0, 5]).
       
Examples:
--------------------------------------------------------------------------------------------------
import Rbeast as rb

Nile, year = rb.load_example('Nile')   # Nile river annual streamflow: trend-only data
metadata   = rb.args()       # an empty object to suff more attribute fields
metadata.season    = 'none'  # trend-only
metadata.startTime = 1871;
metadata.deltaTime = 1;
       
o = rb.beast123(Nile,metadata)   # Default values will be used if parameters are missing
rb.print(o)
rb.plot(o)
 
ohio = rb.load_example('ohio')     # irregular Landsat NDVI time series
metadata = lambda:None;            # equivalent to 'rb.args()', creating an empty dummy object
metadata.season     = 'harmonic'  # seasonality is present
metadata.time       = ohio.time;
metadata.deltaTime  = 1/12;        # aggregate into a monthly ts
metadata.period     = 1.0;         # period= 1/12 x 12 (freq)

o = rb.beast123(ohio.ndvi,metadata)   # Default values will be used if parameters are missing
rb.print(o)
rb.plot(o, ncpStat='median')       
 
ohio.datestr1   # strings of times for ohio.ndvi
metadata              = lambda:None
metadata.time         = rb.args( dateStr=ohio.datestr1, strFmt='????yyyy?mm?dd' )
metadata.deltaTime    = 1/12;   # aggregate into a monthly ts
metadata.period       = 1.0;     # period= 1/12 x 12 (freq)
# Default values will be used if parameters are missing
o =  rb.beast123(ohio.ndvi,metadata,[],[], rb.args(dumpInputData=True)) 
rb.plot(o, ncpStat='median')   
 
ohio.datestr2      # strings of times for ohio.ndvi
metadata      = lambda:None
metadata.time         = lambda:None
metadata.time.dateStr = ohio.datestr2;
metadata.time.strFmt  = 'LC8-yyyydoyxdvi';
metadata.deltaTime    = 1/12;   # aggregate into a monthly ts
metadata.period       = 1.0;    # period= 1/12 x 12 (freq)
o = rb.beast123(ohio.ndvi,metadata) 
rb.plot(o, ncpStat='median')  

#See https://rdrr.io/cran/Rbeast/man/beast123.html for more details
#about the accepted formats of date strings
 
simdata = rb.load_example('simdata')      # a 774x3 matrix: 774 is the time dimesnion      
# simdata is a toy example of 3 time series: I decides to be lazy here bcz the 3
# time series are essentially the same. We just assume they are
# different for illustring the use of beast123 to handle 2d matrix      
 
metadata = lambda:None  
metadata.whichDimIsTime = 1   # 774 is the ts length 
o = rb.beast123(simdata ,metadata) 
rb.print(o,  1)            #print the result for the first ts
rb.plot(o, index=1)        #plot the result for the first ts
rb.print(o,  2)            #print the result for the 2nd ts
rb.plot(o, index=2)        #plot the result for the 2nd ts      
 
import numpy as np
simData3x774 = np.transpose( simdata ) # a 3x774 matrix: 774 is the time dimesnion      
metadata = lambda:None
metadata.whichDimIsTime = 2 # 774 is the ts length 
o = rb.beast123(simData3x774 ,metadata) 
rb.print(o,  2) #print the result for the 2nd ts
rb.plot(o, index=2)     #plot the result for the 2nd ts   
 
ndvi3d, fnames, fyear = rb.load_example( 'imagestack' ) 
# A toy example of stacked time series images: unevely-spaced in time
ndvi3d    # a 12x9x1066 3D cube, and 1066 is the time series length
fnames    # a vector of image filenames as date strings
fyear     # a vector of fractional year 
metadata = lambda:None
metadata.time         = lambda:None ;
metadata.time.dateStr = fnames
metadata.time.strFmt  = 'LT05_018032_20110726.yyyy-mm-dd';
metadata.deltaTime    = 1/12; # aggregated at a monthly interval
metadata.period= 1.0;  # the period is 1.0 (year)
extra   =  rb.args();
extra.dumpInputData    = True # get a copy of the aggregated input
extra.numThreadsPerCPU = 2;  # 2 threads per CPU core
o = rb.beast123(ndvi3d,metadata,[],[], extra) 
rb.print(o,[2,4])      #print the result at row 3 and col 5    
rb.plot(o,index=[2,4]) #plot the result at row 3 and col 5    

metadata = lambda:None
metadata.time         = fyear;
metadata.deltaTime    = 1/12; # aggregated at a monthly interval
metadata.period= 1.0;  # the period is 1.0 (year)
extra =  lambda:None
extra.dumpInputData    = true # get a copy of the aggregated input
extra.numThreadsPerCPU = 2;   # 2 threads per CPU core
o = rb.beast123(ndvi3d,metadata,[],[], extra) 
rb.print(o,[2,4])      #print the result at row  3 and col 5    
rb.plot(o,index=[2,4]) #plot the result at row  3 and col 5     
 
Contact info: To report bug or get help, do not hesitate to contact Kaiguang Zhao
at zhao.1423@osu.edu.
      """
      Y = force_convert_to_numpy(Y)
      
      if hasattr(metadata, 'time'):
            time = metadata.time
            if hasattr(time, "year"):
                time.year = force_convert_to_numpy(time.year)
                if hasattr(time, "month"):
                    time.month = force_convert_to_numpy(time.month)
                if hasattr(time, "day"):
                    time.day   = force_convert_to_numpy(time.day)   
                if hasattr(time, "doy"):
                    time.doy   = force_convert_to_numpy(time.doy)   
            elif hasattr(time,'datestr'):
                time.datestr   = force_convert_to_numpy(time.datestr)  
            elif hasattr(time,'dateStr'):
                time.dateStr   = force_convert_to_numpy(time.dateStr) 
                pass
            else: 
                time = force_convert_to_numpy(time)             #then, we assume time is a numeric vector
                
            metadata.time = time
            
      o = cb.Rbeast('beast_' + method,Y, metadata, prior, mcmc, extra)
      return (o)


