function [Z,W] = computeZW_lambda(marg,Y,A,XX)
Z=[];
W=[];
N=length(A);
for i=1:length(marg)
    W1 = kr(A{marg{i}(2)},A{marg{i}(1)});
    W=[W; W1];
    Z = [Z; Y{i}(:)];    
end
for i=1:N
    W1 = kr(A{i},A{i});
    W=[W; W1];
    Z = [Z; XX{i}(:)];  
end
end