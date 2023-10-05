# Rbeast: A Python package for Bayesian changepoint detection and time series decomposition
 
####  BEAST (Bayesian Estimator of Abrupt change, Seasonality, and Trend) is a fast, generic Bayesian model averaging algorithm to decompose time series or 1D sequential data into individual components, such as abrupt changes, trends, and periodic/seasonal variations, as described in <ins>[Zhao et al. (2019)](https://go.osu.edu/beast2019)</ins>. BEAST is useful for *changepoint detection (e.g., breakpoints, structural breaks, regime shifts, or anomalies), trend analysis, time series decomposition (e.g., trend vs seasonality), time series segmentation, and interrupted time series analysis*. See a list of <a href="#publications"> selected studies using BEAST </a>.

**Quick Installation**
> BEAST was impemented in C/C++ but accessible from  R, Python, Matlab, and Octave.  Run the following to install:

* [Python](#python):   **`pip install Rbeast`**   
* [Matlab](#matlab):  **`eval(webread('http://b.link/rbeast',weboptions('cert','')))`**
* [Octave](#octave):  **`eval(webread('http://b.link/rbeast'))`**  
* [R lang](#r):  **`install.packages("Rbeast")`** 

   
**Quick Usage**

> One-liner code for Python, Matlab and R. Check [github.com/zhaokg/Rbeast](https://github.com/zhaokg/Rbeast) for more details.
```
# Python example
import Rbeast as rb; (Nile, Year)=rb.load_example('nile'); o=rb.beast(Nile,season='none'); rb.plot(o)

# Matlab/Octave example
load('Nile'); o = beast(Nile, 'season','none'); plotbeast(o)

# R example
library(Rbeast); data(Nile); o = beast(Nile); plot(o)
```



## Installation for Python

<p  align="left">   
 <a href= "https://github.com/zhaokg/Rbeast"> <img src="https://img.shields.io/static/v1?style=plastic&logo=github&label=see also&message=github.com/zhaokg/Rbeast&color=brightgreen" height="20"></a>
</p> 

A package **`Rbeast`** has been deposited here at PyPI: https://pypi.org/project/Rbeast/. Install from a binary wheel using:
 
  ```python
    pip install Rbeast
  ```
  Currently, binary wheel files were built for common OS platforms and CPU architectures (e.g., Linux, Windows, and macOS for both x86_64 and arm64 CPUs). If the pre-built wheel doesn't work on your computer, please try to install from the source:
 
  ```python
  pip install Rbeast --no-binary :none:
  ```
 The installation from sources requires a C/C++ compiler (e.g., gcc, and clang) to build the binary package, which should be hassle-free on Linux ( with gcc) or Mac (with clang or xcode) and may be tricky on Windows systems. If needed, contact Kaiguang Zhao (zhao.1423@osu.edu) to help build the package for your specific OS platforms and Python versions.

 ## Run and test Rbeast in Python

Import the Rbeast package as `rb`:
  ```python
import Rbeast as rb
help(rb.load_example)
help(rb.beast)
  ```
The first example is annual streamflow of the River Nile starting from Year 1871. As annual observations, it has no periodic component (i.e., `season='none'`).
  ```python
nile, year = rb.load_example('nile')                     # nile is a 1d Python array or numpy vector
o          = rb.beast( nile, start=1871, season='none')  # season='none' bcz the data has no seasonal/periodic component
rb.plot(o, title='Annual streamflow of the Nile River')
rb.print(o)

# Print a list of fields in the output variable (e.g, o.data, o.RMSE, o.trend.cp, o.time, and o.tend.cpOccPr)
# Check the R manual for expalanations of the output (https://cran.r-project.org/web/packages/Rbeast/Rbeast.pdf) 
o                                                        # this is equivalent to "print(o)"                                    
```

```python
# A wrong way to call beast for nile: If the 'season' arg is missing, the default season='harmonic' is used such that
# there is a seasonal component to be fit. But the Nile data is a trend-only signal with no periodic component
o = rb.beast(nile, start=1871 )  
rb.plot(o)      # the result is wrong. Use season='none' when calling beast for trend-only data                                 
```
  ![](  https://github.com/zhaokg/Rbeast/raw/master/Python/Nile.png)

 The second example is a monthly time series of the Google Search popularity of `beach` over the US. This time series is reguarly-spaced (i.e., deltat=`1 month` =`1/12 year`); it has a cyclyic component with a period of 1 year. That is, the number of data points per period is `period / deltat` =  1 year / 1 month = 1/(1/12) = 12.
 

  ```python
beach, year = rb.load_example('googletrend')
o = rb.beast(beach, start= 2004, deltat=1/12, period = 1.0)       # the time unit is uknown or arbitrary
o = rb.beast(beach, start= 2004, deltat=1/12, period ='1.0 year') # the time unit is fractional year
rb.plot(o)
rb.print(o)
  ```
  ![](  https://github.com/zhaokg/Rbeast/raw/master/Python/beach.png)
  The third example is a stack of 1066 satellite NDVI images over time, with a spatial dimenion of 10 rows x 9 cols: Each pixel is an IRREGULAR time series of 1066 NDVI values with periodic variations at a period of 1.0 year. When running, BEAST will first aggragate the irregular time series into regular ones at a specified time interaval of `deltat` (in this example, we choose `deltat`=`1/12 year` =`1 month`, but you may choose other intervals, depending on the needs).

  ```python 
ndvi3d, datestr,year = rb.load_example('imagestack')  # ndvi is a numpy array of shape (484,10,20); the 1st dim refers to the time

metadata                = rb.args() # create an empty object to stuff the attributes: "metadata  = lambda: None" also works
metadata.time           = year      # times of individulal images/data points: the unit here is fractional year (e.g., 2004.232)
metadata.deltaTime      = 1/12      # regular interval used to aggregate the irregular time series (1/12 = 1/12 year = 1 month)
metadata.period         = 1.0       # the period is 1.0 year, so freq= 1.0 /(1/12) = 12 data points per period
metadata.whichDimIsTime = 1         # the dimension of the input ndvi is (484,10,20): which dim refers to the time. whichDimIsTime is a 1-based index  

o = rb.beast123(ndvi3d, metadata, [], [], []) # rb.beast123(data, metadata, prior, mcmc, extra): default values used if not supplied

rb.print(o[4, 6])                 # print the (5-th row, 7-th col) pixel: Python uses 0-based indices.
rb.plot(o[4, 6])                  # plot the (5-th row, 7-th col) pixel: Python uses 0-based indices.

figure, axes = rb.plot(o[4, 6])   # plot the (5-th row, 7-th col) pixel: Python uses 0-based indices.
rb.plot( o[4, 7], fig = figure)   # plot the (5-th row, 8-th col)pixel: Setting fig=figure will use the existing figure to plot
  ```

Below is another way to supply the time info:
    
  ```python 
ndvi3d, datestr,year = rb.load_example('imagestack') 

metadata              = lambda: None # create an empty object to stuff the attributes: "metadata  = rb.args()" also works
metadata.time         = rb.args( )   # create an empty object to stuff the 'datestr' and 'strfmt' attributes
metadata.time.datestr = datestr      # datestr is a list of file names (e.g., s2_ndvi_2018-01-03.tif) that contain the date info
metadata.time.strfmt  = 'LT05_018032_20110726.yyyy-mm-dd'  # the format used to extract the year (YYYY), month (mm), and day (dd) from the strings
metadata.deltaTime    = 1/12        # regular interval used to aggregate the irregular time series (1/12 = 1/12 year = 1 month)
metadata.period       = 1.0         # the period is 1.0 year, so freq= 1.0 /(1/12) = 12 data points per period
metadata.whichDimIsTime = 1         # the dimension of the input ndvi is (484,10,20): which dim refers to the time. whichDimIsTime is a 1-based index  

extra = rb.args(                             # a set of options to specify the outputs or computational configurations
                 dumpInputData    = True,    # make a copy of the aggregated input data in the beast ouput
                 numThreadsPerCPU = 2,       # Paralell  computing: use 2 threads per cpu core
                 numParThreads    = 0        # `0` means using all CPU cores: total num of ParThreads = numThreadsPerCPU * core Num           
                )
# Instead of using "extra=lambda:None;  extra.dumpInputData=True; ...", the above directly specifies the attribues in the object creation function
 
o = rb.beast123(ndvi3d, metadata, [], [], extra)  # rb.beast123(data, metadata, prior, mcmc, extra): default values used for prior and mcmc if missing
```

## Description
Interpretation of time series data is affected by model choices. Different models can give different or even contradicting estimates of patterns, trends, and mechanisms for the same data–a limitation alleviated by the Bayesian estimator of abrupt change,seasonality, and trend (BEAST) of this package. BEAST seeks to improve time series decomposition by forgoing the "single-best-model" concept and embracing all competing models into the inference via a Bayesian model averaging scheme. It is a flexible tool to uncover abrupt changes (i.e., change-points), cyclic variations (e.g., seasonality), and nonlinear trends in time-series observations. BEAST not just tells when changes occur but also quantifies how likely the detected changes are true. It detects not just piecewise linear trends but also arbitrary nonlinear trends. BEAST is applicable to real-valued time series data of all kinds, be it for remote sensing, finance, public health, economics, climate sciences, ecology, and hydrology. Example applications include its use to identify regime shifts in ecological data, map forest disturbance and land degradation from satellite imagery, detect market trends in economic data, pinpoint anomaly and extreme events in climate data, and unravel system dynamics in biological data. Details on BEAST are reported in [Zhao et al. (2019)](https://go.osu.edu/beast2019). The paper is available at https://go.osu.edu/beast2019.

## Reference
* Zhao, K., Wulder, M. A., Hu, T., Bright, R., Wu, Q., Qin, H., Li, Y., Toman, E., Mallick B., Zhang, X., & Brown, M. (2019). [Detecting change-point, trend, and seasonality in satellite time series data to track abrupt changes and nonlinear dynamics: A Bayesian ensemble algorithm.](https://go.osu.edu/beast2019) Remote Sensing of Environment, 232, 111181. (the BEAST paper) 

* Zhao, K., Valle, D., Popescu, S., Zhang, X. and Mallick, B., 2013. [Hyperspectral remote sensing of plant biochemistry using Bayesian model averaging with variable and band selection](https://www.academia.edu/download/55199778/Hyperspectral-biochemical-Bayesian-chlorophyll-carotenoid-LAI-water-content-foliar-pigment.pdf). Remote Sensing of Environment, 132, pp.102-119. (the mcmc sampler used for BEAST)

* Hu, T., Toman, E.M., Chen, G., Shao, G., Zhou, Y., Li, Y., Zhao, K. and Feng, Y., 2021. [Mapping fine-scale human disturbances in a working landscape with Landsat time series on Google Earth Engine](https://pages.charlotte.edu/gang-chen/wp-content/uploads/sites/184/2021/05/Hu_2021_BEAST-HF-s.pdf). ISPRS Journal of Photogrammetry and Remote Sensing, 176, pp.250-261. (an application paper)


----
<a name=publication></a>

<h2 id="publicationid"> Selected publications using BEAST/Rbeast  </h2> 
 
 Despite being published originally for ecological and enviornmental applications, BEAST is developed as a generic tool applicable to time series or time-series-like data arising from all disciplines. BEAST is not a heuristic algorithm but a rigorous statistical model. Below is a short list of peer-reviewed pulications that used BEAST for statistical data analysis.
 


## Reporting Bugs or getting help

BEAST is distributed as is and without warranty of suitability for application. The one distributed above is still a beta version, with potential room for further improvement. If you encounter flaws/bugs with the software, please report the issue. Providing a description of the conditions under which the bug occurred will help to identify the bug. You can directly email its maintainer Dr. Kaiguang Zhao at <zhao.1423@osu.edu> to report issues or request feature enhancements. Alternatively, use the [Issues tracker](https://github.com/zhaokg/Rbeast/issues) on GitHub. 
