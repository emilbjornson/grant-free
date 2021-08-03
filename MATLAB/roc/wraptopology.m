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


function [dist] = wraptopology(AP,UT,D)
%Here we consider a wrap toplogy and try to find the minimum distance
%between the AP and UT in wrap topology network rather than considering
%only one Cellular Netowork, which would be practical.

M=size(AP,1);
K=size(UT,1);

D1=zeros(M,2);
D1(:,1)=D1(:,1)+ D*ones(M,1);
AP1=AP+D1;

D2=zeros(M,2);
D2(:,2)=D2(:,2)+ D*ones(M,1);
AP2=AP+D2;

D3=zeros(M,2);
D3(:,1)=D3(:,1)- D*ones(M,1);
AP3=AP+D3;

D4=zeros(M,2);
D4(:,2)=D4(:,2)- D*ones(M,1);
AP4=AP+D4;

D5=zeros(M,2);
D5(:,1)=D5(:,1)+ D*ones(M,1);
D5(:,2)=D5(:,2)- D*ones(M,1);
AP5=AP+D5;

D6=zeros(M,2);
D6(:,1)=D6(:,1)- D*ones(M,1);
D6(:,2)=D6(:,2)+ D*ones(M,1);
AP6=AP+D6;

D7=zeros(M,2);
D7=D7+ D*ones(M,2);
AP7=AP+D7;

D8=zeros(M,2);
D8=D8- D*ones(M,2);
AP8=AP+D8;

dist=zeros(M,K);
for m=1:M
    for k=1:K
        dist(m,k) = min([norm(AP(m,:)-UT(k,:)), 
                        norm(AP1(m,:)-UT(k,:)), 
                        norm(AP2(m,:)-UT(k,:)),
                        norm(AP3(m,:)-UT(k,:)),
                        norm(AP4(m,:)-UT(k,:)),
                        norm(AP5(m,:)-UT(k,:)),
                        norm(AP6(m,:)-UT(k,:)),
                        norm(AP7(m,:)-UT(k,:)),
                        norm(AP8(m,:)-UT(k,:)) ]);  
    end
end
