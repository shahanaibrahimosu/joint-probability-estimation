function [Labels] = predict_parafac_multiclass(A,Dataset,I,F)
N      = length(I);
s_test = size(Dataset,1);
Labels = zeros(s_test,1);

num_classes = size(A{1},1);

for n=1:N
    %A{n}(A{n}==0)=1;
    A{n}=5.*A{n};
end
for s = 1:s_test
    p = ones(num_classes,F);
    for n = 2:N
        for n_class = 1:num_classes
            if (Dataset(s,n-1)>0)
                p(n_class,:) =  p(n_class,:).*A{n}(Dataset(s,n-1),:);
            end
        end
    end
    
    for n_class = 1:num_classes
        p(n_class,:) =  p(n_class,:).*A{1}(n_class,:);
    end
    p;
    [~,Labels(s)] = max(sum(p,2));
end
end

