function [A_est1,l_est] = SingleSPA(Y,I,M,F_true,opts)
% input:
%       y   : measurement vector
%       I   : number of rows of each factor
%       opts.
%           constraint  = 'nonnegative'
%                       = 'simplex_col'
%                       = 'simplex'
%           max_iter    : max number of iterations, default: 1000
%           tol_impr    : min cost improvement to terminate the algorithm, default: 1e-7
% Output:
%       A   : cell array containing the factors
%       Out.
%           iter            : number of iterations
%           hist_cost       : history of cost
%           hist_rel_cost   : history of relative cost
%           time_instants   : time at the end of each iteration
%%
%Construct the matrix X to be factorised using the marginal PMF
X_cell =get_golden_matrix(Y,opts.marg,M,I);
X=X_cell{1,1};
X_sum=sum(X,1);
X_n = X*diag(1./sum(X,1));
[l,jk,kl]=FastSepNMF(X_n,F_true,0);
W_est= X(:,l);

Z = X';
W = W_est';
%A = 1/F_true*(ones(size(Z,1),size(W,1)));
A = rand(size(Z,1),size(W,1));
A=diag(1./sum(A,2))*A;
epsilon=1/max(eig((W)*(W')));

for i=1:opts.max_iter
    A = A- epsilon* (A*(W)*(W')-Z*W');
    %A = ProjectOntoSimplex(A',1);
    A=max(A,0);
    %A=A';
    if((norm(Z-A*W)/length(X))<opts.tol)
        break;
    end
end
%H_est=(pinv(W_est)*X);
%H_est = H_est';
%H_est = diag(transpose(X_sum))*H_est;
H_est=A;
N = length(I);
A_est = cell(N,1);
A_est1 = cell(N,1);
W_f = [];
for n=1:M
    A_est{n} = W_est(sum(I(1:n-1))+1:sum(I(1:n)),:);
    A_est1{n} = A_est{n}*diag(1./sum(A_est{n},1));
    W_f = [W_f;A_est1{n}];
end


temp=W_est./W_f;
temp(isnan(temp))=0;
l_est1=max(temp,[],1)';

H_f =[];
for n=M+1:N
    A_est{n} = H_est(sum(I(M+1:n-1))+1:sum(I(M+1:n)),:);
    A_est1{n} = A_est{n}*diag(1./sum(A_est{n},1));
    H_f = [H_f;A_est1{n}];
end

temp=H_est./H_f;
temp(isnan(temp))=0;
l_est2=max(temp,[],1)';  

l_est = l_est1.*l_est2;
l_est = l_est./sum(l_est,1);

end
