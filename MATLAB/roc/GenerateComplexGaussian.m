function ret_val = GenerateComplexGaussian(row,col,var)

    ret_val = sqrt(var)*(1/sqrt(2)*complex(randn(row,col),randn(row,col))) ;

end