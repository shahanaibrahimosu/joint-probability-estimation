function [Z,W] = computeZ(marg,rows,cols,lambda,Y,A)
Z=[];
W=[];

for i=1:length(rows)
    W1 = diag(lambda)*A{marg{rows(i)}([1:cols(i)-1 cols(i)+1:end])};
    W=[W W1];
    if(cols(i)==1)
        Z1=Y{rows(i)};
    else
        Z1 = Y{rows(i)}';
    end
    Z = [Z Z1];    
end
end