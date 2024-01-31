function passed = goodness_of_fit_test(f, x, df, sl)
%f...model function, evaluated for the estimated parameters
%x...measurement (image of molecule)
%sl...significance level; e.g. sl=0.05 for a 5% risk that we kick out data that actually fulfilled our null-hypothesis 
%df...degrees of freedom of the chi2-distribution: number of pixels in the
%molecule minus number of fitting parameters, e.g. 15x15 - 5 = 220. 

LLR = 2*sum(sum((f - x .* log(f)) - (x - x .* log(x))));

if chi2cdf(LLR,df,'upper') <= sl %if the probability that a molecule image fulfilling the null hypothesis is more "deviating" than the observed one is smaller than sl
    passed = 0; %fit fails
else 
    passed = 1;
end

end




