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

function [gamma_hat] = decode_activity_pattern_poly(sigma_sqr,Beta,Y,S)
% Function  to decode the activity pattern from two APs

    [K , ~] = size(Beta) ; %  K - #Users 
    [L, N, M] = size(Y) ;  %  L - Sequence length, N - #Antennas per AP, M - #APs
        
    % Initialize the coordinate descent method 
    
    % Initialize Output Vector
    gamma_hat = zeros(K,1) ;
    gamma_prev = gamma_hat ; 
    
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
    
    
    cost_prev = inf ; 
    
    % Iteration index
    for i=1:10
      
        Kidx = randperm(K,K) ; 
        
        % Run through all indices randomly
        for k = Kidx
            m1 = su_link(end,k) ;   % Largest 
            m2 = su_link(end-1,k) ; % Second Largest 
            
            s1_m1 = Qm_inv(:,:,m1)*S(:,k) ;
            s2_m1 = (Qm_inv(:,:,m1)'*S(:,k))' ;
            
            s1_m2 = Qm_inv(:,:,m2)*S(:,k) ;
            s2_m2 = (Qm_inv(:,:,m2)'*S(:,k))' ;
            
            
            N1_m1 = real(s2_m1*QY(:,:,m1)*s1_m1) ; 
            N2_m1 = real(s2_m1*S(:,k)) ; 

            N1_m2 = real(s2_m2*QY(:,:,m2)*s1_m2) ; 
            N2_m2 = real(s2_m2*S(:,k)) ;
            
            
            am1 = Beta(k,m1)*N2_m1 ; 
            am2 = Beta(k,m2)*N2_m2 ; 
            
            bm1 = Beta(k,m1)*N1_m1 ; 
            bm2 = Beta(k,m2)*N1_m2 ; 
            
            % Form the coefficients of Polynomial of degree 3 from am, bm
            c = zeros(1,4) ; 
            c(1) = 2*(am1^2)*(am2^2) ; 
            c(2) = 3*(am1^2)*am2 + 3*(am2^2)*am1 - (am1^2)*bm2 - (am2^2)*bm1 ; 
            c(3) = am1^2 + am2^2 + 2*am1*(am2-bm2) + 2*am2*(am1-bm1) ; 
            c(4) = am1-bm1 + am2-bm2 ; 
            
            % Compute Roots
            r = roots(c) ; 
            
            % Pick Real Valued Roots and add boundary Condition 0
            Val = [r(imag(r)==0)' -1*gamma_hat(k)];  % Added the boundary conditions  
            
            % Look for the minimizer
            if (length(Val) >1)
                f1 = log(1+am1*Val) ; 
                f2 = log(1+am2*Val) ;
                f3 = Val*bm1./(1+am1*Val) ;
                f4 = Val*bm2./(1+am2*Val) ;
                f = f1+f2-f3-f4 ; 
                [~,idx] = min(f) ; 
                Val = Val(idx) ;
            end
            
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