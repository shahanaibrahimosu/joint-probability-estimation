function [Ys] = get_selected_Y(Y,marg,I)
marg_array =cell2mat(marg);
Ys  = cell(size(marg,1),1);
for m=1:size(marg,1)
    j=marg_array(m,1);
    k=marg_array(m,2);
    l=marg_array(m,3);
    %Get the actual tensor
    X_mat = Y{j,k,l};
    Ys{m}= zeros(I(j),I(k),I(l));
    for t=1:I(j)
      Ys{m}(t,:,:)=reshape(X_mat(:,t),[I(k),I(l)]);
    end
end
end