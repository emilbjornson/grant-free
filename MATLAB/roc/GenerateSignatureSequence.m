function S = GenerateSignatureSequence(L,K)
  
    % L : Signature Length
    % K : Number of Sequences
    % each sequence has unit power 
    S = GenerateComplexGaussian(L,K,1) ; 
   
end