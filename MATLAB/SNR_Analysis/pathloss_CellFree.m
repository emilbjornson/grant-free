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



function [BETAA] = pathloss_CellFree(dist)

M=size(dist,1);
K=size(dist,2);

BETAA = zeros(M,K);

for m=1:M
    for k=1:K
        betadB = -30.5 - 36.7*log10(dist(m,k)*1000) + 4*randn(1,1); 
        
        BETAA(m,k)=10^(betadB/10);
        
    end
end

end