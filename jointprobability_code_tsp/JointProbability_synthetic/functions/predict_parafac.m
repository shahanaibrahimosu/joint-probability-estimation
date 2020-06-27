function [Labels,ct] = predict_parafac(A,lambda,Dataset,I,F)
N      = length(I);
s_test = size(Dataset,1);
Labels = zeros(s_test,1);
ct     = 0;
for s = 1:s_test
    
    p1 = ones(1,F);
    p2 = ones(1,F);
    
    for n = 2:N
        p1 = p1.*A{n}(Dataset(s,n-1),:);
        p2 = p2.*A{n}(Dataset(s,n-1),:);
    end
    
    p1 = p1.*A{1}(1,:);
    p2 = p2.*A{1}(2,:);
    
    Labels(s) = (p2*lambda>p1*lambda) + 1;
    
    % Catch weird cases
    if (p2*lambda==0) && (p1*lambda==0)
        Labels(s) = -1;
        ct = ct+1;
    end
end
end

