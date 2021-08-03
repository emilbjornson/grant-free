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

function [gamma_hat] = decode_activity_pattern_Mpoly(sigma_sqr,Beta,Y,S,NAP)
% Function  to decode the activity pattern from two APs

    [K , ~] = size(Beta) ; %  K - #Users 
    [L, N, M] = size(Y) ;  %  L - Sequence length, N - #Antennas per AP, M - #APs
        
    % Initialize the coordinate descent method 
    
    % Initialize Output Vector
    gamma_hat = zeros(K,1) ; 
    
    % Initialize the cost check 
    gamma_prev = gamma_hat ; 
    cost_prev = inf ; 
    
    % Initialize Covariance of Y 
    QY = zeros(L,L,M) ; 
    for m=1:M
        QY(:,:,m) = (1/N)*(Y(:,:,m)*Y(:,:,m)') ;
    end

    % Initialize the Covariance of Noise - Inverse 
    Qm_inv = zeros(L,L,M) ; 
    for m=1:M
        Qm_inv(:,:,m)   = (1/sigma_sqr)*diag(ones(1,L))  ;
    end
    
    % Find the strongest link for each user
    % value range from [1 M]
    [~, su_link] = sort(Beta.') ; 
    
    % Iteration index
    for i=1:10
      
        Kidx = randperm(K,K) ; 
        
        % Run through all indices randomly
        for k = Kidx
                        
            am = zeros(1,NAP) ;
            bm = zeros(1,NAP) ; 
            for n = 0:1:NAP-1
                m = su_link(end-n,k) ; 
                s1 = Qm_inv(:,:,m)*S(:,k) ;
                s2 = (Qm_inv(:,:,m)'*S(:,k))' ;
            
                N1 = real(s2*QY(:,:,m)*s1) ; 
                N2 = real(s2*S(:,k)) ; 
                
                am(n+1) = Beta(k,m)*N2 ;
                bm(n+1) = Beta(k,m)*N1 ;
            end
            
            c = compute_poly_coefficients(am,bm) ; 
            
            % Compute Roots
            r = roots(c) ; 
            
            % Pick Real Valued Roots and add boundary Condition 0
            Val = [r(imag(r)==0)' -1*gamma_hat(k)];  % Added the boundary conditions  
            
            % Look for the minimizer
            f = zeros(1,length(Val)) ; 
            for n=1:1:NAP
                f1 = log(1+am(n)*Val) ;     
                f2 = Val*bm(n)./(1+am(n)*Val) ;
                f = f+f1-f2 ; 
            end
            [~,idx] = min(f) ;
            
                       
            Val = Val(idx) ;
            
            % Compute delta to Update 
            delta = max(Val,-1*gamma_hat(k)) ;

            % Update gamma_hat with delta
            gamma_hat(k) = gamma_hat(k) + delta ;

            % Update every Covariance matrix with delta 
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
        
        % Compute Cost Function for Convergence Check
        cost = 0 ; 
        for m=1:1:M
            cost = cost + log(det(inv(Qm_inv(:,:,m)))) + trace(Qm_inv(:,:,m)*QY(:,:,m)) ; 
        end
        
        if(cost<cost_prev)
            cost_prev = cost ; 
            gamma_prev = gamma_hat ; 
        else
            gamma_hat = gamma_prev ; 
            break ; 
        end
        
    end

end