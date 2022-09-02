function [bf,p] = ttest_bf(t,n)
%%%% This function 


nu = n-1;
r = .707;

denom = @(g) (1+n.*g.*r.^2).^(-.5) .* (1 + (t.^2) ./ ((1+n.*g.*r.^2)*nu) ).^(-.5*n) .* (2*pi).^(-.5) .* g.^(-3/2) .* exp(-1./(2*g));

nominator  = (1 + (t^2 / nu))^(-.5*n);
denominator = integral(denom,0,Inf);

bf = nominator / denominator;

null_likelihood = bf / (1+bf);

p = 1 - null_likelihood;


end