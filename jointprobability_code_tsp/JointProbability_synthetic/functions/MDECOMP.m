function [A,lambda,Out] = MDECOMP(Y,I,F,opts,A_t,l_t)

N= length(I);
M=3;
marg_array =cell2mat(opts.marg);
Z = get_golden_matrix(Y,opts.marg,M,I);
[U,S,V] = svd(Z{1},'econ');
U_est = U(:,1:F);
for i=1:M
    A_est{i} = U_est(I(i)*(i-1)+1:I(i)*i,:);
end
X=cell(M,1);
y=[];
B=[];
for i=1:M
    j=i+1;
    X{i} = Y{marg_array(:,1)==i & marg_array(:,2)==j};
    y = [y;X{i}*ones(I(j),1)];
    B = [B;kron(ones(I(1),1)',A_est{i})];
end
h= B\y;



end