function [A,l] = gen_PMF_factors(I,F)
% Generate factors
N     = length(I);
A     = cell(N,1);
for n = 1:N
    A{n}     = rand(I(n),F);
    A{n}     = bsxfun(@times,A{n},1./sum(A{n},1));
end

l = rand(F,1); l = l/sum(l);
end