function [Labels,ct] = predict_labels(A,Dataset,I,F)
N      = length(I);
s_test = size(Dataset,1);
Labels = zeros(s_test,1);
ct     = [];
for n=1:N
    %A{n}(A{n}==0)=1;
    A{n}=5.*A{n};
end
for s = 1:s_test

    p1 = ones(1,F);
    p2 = ones(1,F);
    
    for n = 2:N
        p1 = p1.*A{n}(Dataset(s,n-1),:);
        p2 = p2.*A{n}(Dataset(s,n-1),:);
    end
    
    p1 = p1.*A{1}(1,:);
    p2 = p2.*A{1}(2,:);
    
    Labels(s) = (sum(p2)>sum(p1)) + 1;
    
    % Catch weird cases
    if (sum(p2)==0) && (sum(p1)==0)
        Labels(s) = -1;
    else
        ct = [ct s];
    end
end
end

