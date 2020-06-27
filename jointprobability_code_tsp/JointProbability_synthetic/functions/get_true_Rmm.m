function [Y] = get_true_Rmm(A,I,lambda)
Y=cell(size(A,1),1);
for m=1:size(A,1)
    %Get the actual tensor
    X_mat = A{m}*diag(lambda)*A{m}';
    Y{m,1}=X_mat;   
end
end

