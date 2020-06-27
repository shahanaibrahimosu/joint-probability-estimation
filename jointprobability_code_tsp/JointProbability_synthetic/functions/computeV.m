function V = computeV(marg,rows,cols,lambda,A,Y,I)
F = length(lambda); 
V = zeros(I,F);

for i=1:length(rows)
    V = V + mtkrprod(Y{i},A(marg{rows(i)}),cols(i));
end
V = diag(lambda)*V';
