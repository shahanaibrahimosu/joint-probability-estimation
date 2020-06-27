function [Y] = get_true_marg_fast(A,marg,I,lambda)
marg_array =cell2mat(marg);
Y=cell(size(marg,1),1);
for m=1:size(marg,1)
    j=marg_array(m,1);
    k=marg_array(m,2);
    
    %Get the actual tensor
    X_mat = A{j}*diag(lambda)*A{k}';
    Y{m,1}=X_mat;   
end
end

