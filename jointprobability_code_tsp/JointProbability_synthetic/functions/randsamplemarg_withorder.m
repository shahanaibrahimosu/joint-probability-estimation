function [X,ind,k,l,row1] = randsamplemarg_withorder(j,Y,I,marg,order_j)
  marg_array=cell2mat(marg);
  %Sample k and l
  N = length(I);
  N_array = 1:N;
  k_array = N_array(N_array ~= j);
  k = datasample(k_array,1);
  l_array = k_array(k_array ~= k);
  l = datasample(l_array,1);
  
  %Get the corresponding marginal for j,k,l 
  [ind,~]=sort([j,k,l]);
  [row,col] =find(marg_array(:,1)==ind(1) & marg_array(:,2)==ind(2) & marg_array(:,3)==ind(3));
  X = Y{row,col};
  [row1]=find(order_j==row);
end