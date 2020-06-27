function [Y] = get_true_marg_fast(A,marg,I,lambda)
marg_array =cell2mat(marg);
Y=cell(size(marg,1));
for m=1:size(marg,1)
    j=marg_array(m,1);
    k=marg_array(m,2);
    l=marg_array(m,3);
    
    %Get the actual tensor
    X_mat = kr(A{l},A{k})*diag(lambda)*A{j}';
    X_tens= zeros(I(j),I(k),I(l));
    for t=1:I(j)
      X_tens(t,:,:)=reshape(X_mat(:,t),[I(k),I(l)]);
    end
    Y{m}=X_tens;   
end
end

