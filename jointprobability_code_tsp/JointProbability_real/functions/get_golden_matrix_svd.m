function [X] = get_golden_matrix_svd(Y,marg,M,I)
marg_array =cell2mat(marg);
N=length(I);

L=1;
X= cell(L,2);
N_list = 1:N;
for k=1:L
X1=cell(N,1);
for i=N_list(1:N)
    X1{i}=[];
    for j=N_list(1:N)
        if(j<i)
            X1{i} = [X1{i} ; Y{marg_array(:,1)==j & marg_array(:,2)==i}];
        else
            X1{i} = [X1{i} ; Y{marg_array(:,1)==i & marg_array(:,2)==j}'];
        end
    end
    X{k,1}=[X{k,1} X1{i}];
end
%X{k,1} = X{k,1}*diag(1./sum(X{k,1},1));
X{k,2}= N_list(1);
N_list=circshift(N_list,-M);
end
end