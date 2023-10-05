# Rbeast: A Python package for Bayesian changepoint detection and time series decomposition
 
####  BEAST (Bayesian Estimator of Abrupt change, Seasonality, and Trend) is a fast, generic Bayesian model averaging algorithm to decompose time series or 1D sequential data into individual components, such as abrupt changes, trends, and periodic/seasonal variations, as described in <ins>[Zhao et al. (2019)](https://go.osu.edu/beast2019)</ins>. BEAST is useful for *changepoint detection (e.g., breakpoints, structural breaks, regime shifts, or anomalies), trend analysis, time series decomposition (e.g., trend vs seasonality), time series segmentation, and interrupted time series analysis*. See a list of <a href="#publications"> selected studies using BEAST </a>.

**Quick Installation**
> BEAST was impemented in C/C++ but accessible from  R, Python, Matlab, and Octave.  Run the following to install:

* [Python](#abcc):   **`pip install Rbeast`**   
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

<h2 name="publications" id ="publications">aaaa</h2>

<p  align="left">   
 <a href= "https://github.com/zhaokg/Rbeast"> <img src="https://img.shields.io/static/v1?style=plastic&logo=github&label=see also&message=github.com/zhaokg/Rbeast&color=brightgreen" height="20"></a>
</p> 

 

 
