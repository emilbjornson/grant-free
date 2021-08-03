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

function [c] = compute_poly_coefficients(am,bm)
    M = length(am) ; 

    c = zeros(1,2*M) ;
    for m=1:1:M
        temp = 1 ; 
        for k=1:1:M
            if k==m
                continue;
            end
            temp = conv(temp,[am(k)^2 2*am(k) 1]) ;
        end
        temp = conv(temp,[am(m)^2 am(m)-bm(m)]) ;
        c = c +  temp ; 
    end
   
end