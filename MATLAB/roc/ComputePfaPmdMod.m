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



function [Pmd, Pfa] = ComputePfaPmdMod(gamma,gamma_hat,threshold)

    [K, monte] = size(gamma) ; 
    
    Pd = zeros(1,monte);
    Pfa = zeros(1,monte);
    for m=1:1:monte
            A = sum(gamma(:,m)>0) ; % Find number of active elements in each MC simulation
            Pd(m)  = sum(( (gamma(:,m)>0) & (gamma_hat(:,m)>threshold) ))/A  ;
            Pfa(m) = sum( ((gamma_hat(:,m)>threshold) & (~(gamma(:,m)>0)) ))/(K-A)  ;
    end

        Pd  = mean(Pd)  ; 
        Pfa = mean(Pfa) ; 
        Pmd = 1-Pd ;
end