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

% Define Number of APs
M = 20; 

% Number of Antennas at each AP
N = 2   ;  

% Total Number of Users in Cell Network
K = 400 ; 

% Binomial Activity percentage : Activity pattern
epsilon = 0.1 ; 

% Pilot Sequence Length
L = 40  ;  % Pilot Length

% Cell Area, Shadow Fading parameter
D=2; % area size in kilometer
SNR = 6 ; % 6dB for D=2, 17.2dB for D=1 % 28dB for D=0.5

TxPow = 200e-3 ; % Transmit Power Constraint in W

% Noise Parameter : bandwidth = 1 MHz
sigma_sqr_dBm = -109 ;
sigma_sqr = 10^((sigma_sqr_dBm-30)/10) ; 
sigma_sqrN = 1 ; % Normalized Sigma2

% Generate Signature Sequence
S = GenerateSignatureSequence(L,K) ;

monte = 1e5 ; 

% get the variables for PFA and PMD computation
gamma=zeros(K,monte) ; 
gamma_hat=zeros(K,monte) ; 


    parfor i = 1:1:monte

        rho       = zeros(K,1);
                
        AP=unifrnd(-D/2,D/2,M,2); % Uniformly Allocate the AP locations
        UT=unifrnd(-D/2,D/2,K,2); % Uniformly Allocate the UT locations

        dist = wraptopology(AP,UT,D) ; 

        Beta = pathloss_CellFree(dist) ;  % MxK Matrix 

        maxBeta = max(Beta) ; % Find the maximum of beta to have power allocation.
        Ac_list = binornd(1,epsilon,[K,1]) ;
        idx = find(Ac_list) ;
        
        SNR_BS = 10*log10(TxPow*maxBeta(idx)/sigma_sqr) ; 
        
        % Drop the users who doesnt meet SNR constraints
        idx = idx(SNR_BS>=SNR) ; 
        
        Beta = Beta' ; % make the size KxM 

        % Power Control for best AP. 
        % Cut the power such that SNR is same for all devices at the best
        % AP
        rho(idx) =  10^(SNR/10)./maxBeta(idx) ;

        gamma(:,i) = rho.*maxBeta' ;
 
        
        X = diag(sqrt(rho)) ;

        Y = zeros(L,N,M) ;

        % Find the signal receieved by each APs
        for m=1:M
            H = (1/sqrt(2)*complex(randn(K,N),randn(K,N))) ;
            W = (1/sqrt(2)*complex(randn(L,N),randn(L,N))) ;
            G = sqrt(Beta(:,m)).*H ; 
            Y(:,:,m) = S*X*G + W ; 
        end

%         Decode with 1 APs        
%         gamma_hat(:,i) = su_decode_activity_pattern(sigma_sqrN,Beta,Y,S).*maxBeta' ; 

%         Decode with 2 APs
        gamma_hat(:,i) = decode_activity_pattern_poly(sigma_sqrN,Beta,Y,S).*maxBeta' ; 

%         Decode with 3 AP
%         gamma_hat(:,i) = decode_activity_pattern_Mpoly(sigma_sqrN,Beta,Y,S,3).*maxBeta' ; 

      
%         Algorithm 3 Parallel Computation        
%         gamma_hat(:,i) = decode_activity_pattern_Mpoly_LC(sigma_sqrN,Beta,Y,S,2,400).*maxBeta' ; 

    end
    
close all
    
thrMin = 0 ; 
thrMax = max(max(gamma)) ; 

threshold = linspace(thrMin,thrMax,200) ; 

PFA = zeros(1,length(threshold)) ;
PMD = zeros(1,length(threshold)) ;
   
for i=1:length(threshold)
    [Pmd,Pfa] = ComputePfaPmdMod(gamma,gamma_hat,threshold(i)) ;
    PMD(i) = Pmd;
    PFA(i) = Pfa;
end


PFA_PMD = [PFA' PMD'] ; 
%save('CellFree_CD_PFA_PMD_M20_N2_K400_L40_D2_200mW_2AP.mat','PFA_PMD');

loglog(PFA,PMD)
grid on ;

set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultTextInterpreter','latex')


ylabel('Probability of Miss Detection');
xlabel('Probability of False Alarm');


set(gca,'FontSize',13)
%savefig('CellFree_CD_PFA_PMD_M20_N2_K400_L40_D2_200mW_2AP.fig')
