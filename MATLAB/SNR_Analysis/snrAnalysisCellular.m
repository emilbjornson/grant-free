%This Matlab function was developed to generate simulation results to:
%
%Unnikrishnan Kunnath Ganesan, Emil Björnson and Erik G. Larsson (2021), 
%"Clustering Based Activity Detection Algorithms for Grant-Free Random
% Access in Cell-Free Massive MIMO", 
%IEEE Transactions in Communications
%
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.


function [SNR] = snrAnalysisCellular(D,TxPow,sigma_sqr,monte_samples)

SNR = zeros(1,monte_samples) ; 

parfor monte = 1:1:monte_samples

    % Consider AP at the middle of the cell
    AP=[0,0] ; 
    UT=unifrnd(-D/2,D/2,1,2); % Uniformly Allocate the UT locations

    dist = wraptopology(AP,UT,D) ; 

    Beta = pathloss_CellFree(dist) ; 

    SNR(monte) = 10*log10(TxPow*Beta/sigma_sqr) ; 
    
end

end