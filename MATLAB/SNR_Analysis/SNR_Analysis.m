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



clc
clear all
close all
set(gca,'DefaultLineLineWidth',2)


% Define Number of APs
M = 20; 

% Number of Antennas at each AP
N = 2   ;  

% Total Number of Users in Cell Network
K = 400 ; 

% Cell Area, Shadow Fading parameter
D=1; % area size in kilometer
    
% Compute Noise Power
sigma_sqr_dBm = -109; % For 1MHz bandwidth
sigma_sqr = 10^((sigma_sqr_dBm-30)/10) ; 

% Transmit Power Limitation in Watts
TxPow = 200e-3 ; 

monte_samples = 1e5 ; 

D = 2 ; 
SNR = snrAnalysisCellular(D,TxPow,sigma_sqr,monte_samples); 
cdfplot(SNR)
hold on

D = 1 ; 
SNR = snrAnalysisCellular(D,TxPow,sigma_sqr,monte_samples); 
cdfplot(SNR)
hold on

D = 2 ; 
SNR = SnrAnalysisCellFree(D,M,TxPow,sigma_sqr,monte_samples); 
cdfplot(SNR)
hold on

D = 1 ; 
SNR = SnrAnalysisCellFree(D,M,TxPow,sigma_sqr,monte_samples); 
cdfplot(SNR)
hold on


grid on
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultTextInterpreter','latex')

xlabel('SNR$_k$ [dB]','Interpreter','latex');
ylabel('CDF','Interpreter','latex');


set(gca,'FontSize',20)
set(gcf,'Position',[200 200 800 600])

xlim([-20 80])

set(findall(gca, 'Type', 'Line'),'LineWidth',1.5);