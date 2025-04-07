## A Parameter and Flag Adaptive Reconstruction Method (PF-Free) for Satellite Vegetation Index Time Series
## Code Author
Yuxi Ran (Wuhan University, Wuhan, China)
## Email
yx.ran@whu.edu.cn

## Brief Introduction of PF-Free
##
The data quality issue induced by atmospheric and other disturbances can significantly impede the application of remote sensing vegetation indices (VIs). 
Despite the development of numerous VI reconstruction techniques, two major challenges remain, i.e., the reliance on quality flag inputs and parameter 
settings. Quality flags are usually necessary as inputs for the different methods to improve the reconstruction accuracy, but mislabeling can be common 
in the quality flag data, which can directly introduce uncertainties. Furthermore, constant parameter schemes are usually assigned during the 
reconstruction applications, but the optimal parameters for all the models can show great spatial heterogeneity. Accordingly, in this paper, we propose
a parameter-free and flag-free adaptive time-series method based on a variational reconstruction framework (PF-Free) to address the afore-mentioned 
issues, which can be applied without any parameter or flag inputs. PF-Free makes full use of the time-series temporal smoothness and inter-annual 
similarity to label the data after time-series rearrangement, and the parameters are adaptively selected using an improved generalized cross-validation 
(GCV) technique. Simulation and real-data experiments all demonstrate that PF-Free can achieve better and more stable reconstruction results, compared
to other comparative methods. The adaptive quality flags can denote the data quality robustly and accurately, while guaranteeing better reconstruction
performance, as long as there is any mislabeling in the original flags. Moreover, the adaptive parameter selection strategy considers the great spatial
heterogeneity in the optimal parameters on a pixel-by-pixel basis, leading to more stable reconstruction outcomes under complex conditions. Further 
experiments also prove the effectiveness of PF-Free in processing data without quality flags or severely contaminated data, using Advanced Very 
High-Resolution Radiometer (AVHRR) data and Moderate Resolution Imaging Spectroradiometer (MODIS) daily normalized difference vegetation index (NDVI) 
data. This work provides a practical reconstruction method for VI time series, which is both flexible and convenient, without requiring any parameter
or flag input, which we believe will advance VI reconstruction applications.
##

## Code Language
MATLAB R2022a

## How to use PF-Free
We have provided a test case (see test.m) with clear annotations

## For more details,see the following papers:
Shen H, Ran Y, Guan X, et al. A Parameter and Flag Adaptive Reconstruction Method (PF-Free) for Satellite Vegetation Index Time Series. lEEE Transactions on Geoscience and Remote sensing.
Craven P, Wahba G. Smoothing noisy data with spline functions: estimating the correct degree of smoothing by the method of generalized cross-validation. Numerische mathematik, 1978, 31(4): 377-403.
Brezinski C, Redivo-Zaglia M, Rodriguez G, et al. Multi-parameter regularization techniques for ill-conditioned linear systems. Numerische Mathematik, 2003, 94: 203-228.
储栋,管小彬,沈焕锋.考虑全时间序列信息的NDVI变分重建方法.地理与地理信息科学,2023,39(03):31-39.
