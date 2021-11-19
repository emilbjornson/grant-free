Clustering-Based Activity Detection Algorithms for Grant-Free Random Access in Cell-Free Massive MIMO
=====================================================================================================

This is a code package is related to the following scientific article:

U. K. Ganesan, E. Björnson and E. G. Larsson, "Clustering-Based Activity Detection Algorithms for Grant-Free Random Access in Cell-Free Massive MIMO," 
in IEEE Transactions on Communications, vol. 69, no. 11, pp. 7520-7530, Nov. 2021, doi: 10.1109/TCOMM.2021.3102635.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results in the article. We encourage you to also perform reproducible research!


## Abstract of Article
Future wireless networks need to support massive machine type communication (mMTC) where a massive number of devices accesses the network and massive MIMO is a promising enabling technology. 
Massive access schemes have been studied for co-located massive MIMO arrays. 
In this paper, we investigate the activity detection in grant-free random access for mMTC in cell-free massive MIMO networks using distributed arrays. 
Each active device transmits a non-orthogonal pilot sequence to the access points (APs) and the APs send the received signals to a central processing unit (CPU) for joint activity detection. 
The maximum likelihood device activity detection problem is formulated and algorithms for activity detection in cell-free massive MIMO are provided to solve it. 
The simulation results show that the macro diversity gain provided by the cell-free architecture improves the activity detection performance compared to co-located architecture when the coverage area is large.

## Content of Code Package

The code package contains codes to plot ROC curves generated in the paper. Use main.m to generate the figures. 
The variables can be configured to get different plots.

See each file for further documentation.

## Acknowledgements

Unnikrishnan Kunnath Ganesan and Erik G. Larsson were supported in part by ELLIIT and in part by Swedish Research Council(VR).
Emil Björnson was supported by the Grant 2019-05068 from the Swedish Research Council.

## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
