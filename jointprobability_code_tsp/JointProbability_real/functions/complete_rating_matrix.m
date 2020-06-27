function [my_rmse,my_mae,pr_rating,num_nan] = complete_rating_matrix(Dataset_mis_val,P4,r4,c4,A,lambda)
pr_rating = zeros(length(P4),1);
nann = 0;

count   = 0;
my_rmse = 0;
my_mae  = 0;

for i=1:length(P4)
    ind_nz    = find(Dataset_mis_val(r4(i),:));
    tmp = ones(1,length(lambda));
    
    if ~isempty(ind_nz)
        for ii = 1:length(ind_nz)
            tmp = tmp.* A{ind_nz(ii)}(Dataset_mis_val(r4(i),ind_nz(ii)),:);
        end   
        tmp    = tmp.*lambda;
        subPMF = A{c4(i)}*tmp';
        subPMF1=sum(subPMF,2);
        subPMF1=subPMF1/sum(subPMF1);
        pr_rating(i) = [1 2 3 4 5]*subPMF1;%*subPMF(:)/sum(subPMF(:));
        
        if isnan(pr_rating(i))
            nann = nann+1;
            pr_rating(i) = 5.5/2;
        end
        
        my_rmse = my_rmse + (pr_rating(i) - P4(i))^2;
        my_mae  = my_mae + abs(pr_rating(i) - P4(i));
        count   = count + 1;
    end
end

fprintf('Number of ratings predicted  : %d  \n', count);
my_rmse = sqrt(my_rmse/count);
my_mae  = my_mae/count;
num_nan = nann;