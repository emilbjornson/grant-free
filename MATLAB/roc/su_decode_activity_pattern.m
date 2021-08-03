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


function [gamma_hat] = su_decode_activity_pattern(sigma_sqr,Beta,Y,S)
% Funtion to decode the activity pattern from the strongest user 

    [K , ~] = size(Beta) ; 
    [L, N, M] = size(Y) ; 
    
    
    % Initialize the coordinate descent method 
    
    gamma_hat = zeros(K,1) ; 
    
    QY = zeros(L,L,M) ; 
    for m=1:M
        QY(:,:,m) = (1/N)*(Y(:,:,m)*Y(:,:,m)') ;
    end

    Qm_inv = zeros(L,L,M) ; 
    for m=1:M
        Qm_inv(:,:,m)   = (1/sigma_sqr)*diag(ones(1,L))  ;
    end
    
    % Find the strongest link for each user
    % value range from [1 M]
    [~, su_link] = max(Beta.') ; 
    
    
    % Iteration index
    for i=1:10
      
        Kidx = randperm(K,K) ; 
        
        % Run through all indices randomly
        for k = Kidx
            m = su_link(k) ; 
            
            s1 = Qm_inv(:,:,m)*S(:,k) ;
            s2 = (Qm_inv(:,:,m)'*S(:,k))' ;
            
            N1 = real(s2*QY(:,:,m)*s1) ; 
            N2 = real(s2*S(:,k)) ; 

            Val = (N1-N2)/(N2^2) ;
            Val = Val/Beta(k,m) ; 
            

            delta = max(Val,-1*gamma_hat(k)) ;

            gamma_hat(k) = gamma_hat(k) + delta ;


            % Update every Covariance matrix with delta from strongest user
            for m=1:M 
                s1 = Qm_inv(:,:,m)*S(:,k) ;
                s2 = (Qm_inv(:,:,m)'*S(:,k))' ;
            
                N2 = real(s2*S(:,k)) ; 
                
                N3 = s1*s2;
                N4 = 1 + (delta*Beta(k,m)*N2);
                q = (delta*Beta(k,m))/N4;
            
                Qm_inv(:,:,m) = Qm_inv(:,:,m) - q*N3 ; 
            end
            
        end
        
    end

end