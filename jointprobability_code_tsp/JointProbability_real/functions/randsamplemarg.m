function [X,ind,k] = randsamplemarg(j,Y,I,marg)
  marg_array=cell2mat(marg);
  %Sample k and l
  N = length(I);
  N_array = 1:N;
  k_array = N_array(N_array ~= j);
  k = datasample(k_array,1);
 
  
  %Get the corresponding marginal for j,k,l 
  [ind,~]=sort([j,k]);
  X = Y{marg_array(:,1)==ind(1) & marg_array(:,2)==ind(2)};
  if(ind(1)~=j)
      X=X';
  end
end