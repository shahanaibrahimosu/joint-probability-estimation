function [A,lambda,Out] = N_CTF_AO_KL(Y,I,F,opts)
% input:
%       y   : measurement vector
%       I   : number of rows of each factor
%       F   : rank
%       opts.
%           max_iter    : max number of iterations, default: 1000
%           tol_impr    : min cost improvement to terminate the algorithm, default: 1e-7
% Output:
%       A   : cell array containing the factors
%       Out.
%% Parameters
N  = length(I);

if ~isfield(opts,'max_iter'),   opts.max_iter   = 1000;          end
if ~isfield(opts,'tol_impr');   opts.tol_impr   = 1e-7;          end

A      = opts.A0;     % Initial tensor factors
lambda = opts.l0;     % Prior of hidden variable
tol    = opts.tol_impr;

% cost and relative cost
Out.hist_cost         = zeros(opts.max_iter,1);
Out.time_instants     = zeros(opts.max_iter,1);

iter = 1;
rows = cell(N,1);
cols = cell(N,1);

for n = 1:N
    for i=1: size(opts.marg,1)
        [~,c]= find(opts.marg{i} == n);
        %if(~isempty(c) && sum(Y{i},"all")> 0)
        if(~isempty(c) && rank(Y{i})> 0)
            rows{n}= [rows{n} i];
            cols{n}= [cols{n} c];
        end
    end
end
Y_n={};
marg_p={};
count=0;
for i=1: size(opts.marg,1)
    if(rank(Y{i})> 0)
        count=count+1;
        Y_n{count,1}=Y{i};
        marg_p{count,1}=opts.marg{i};
    end
end



%%%%%%%% Precompute Data %%%%%%%%
Y_p = cell(N,1);
for n = 1:N
    Y_p{n} = cell(length(rows{n}),1);
    for i=1:length(rows{n})
        Y_p{n}{i} = reshape(permute(Y{rows{n}(i)},[cols{n}(i) 1:cols{n}(i)-1 cols{n}(i)+1:ndims(Y{rows{n}(i)})]),I(n),[]);
        %Y_p{n}{i}(Y_p{n}{i}==0)=1e-6;
        %Y_p{n}{i}=Y_p{n}{i}./sum(Y_p{n}{i},"all");
    end
end
n_tic = tic;
while(1)
    %%%  Solve each subproblem with Mirror Descent
    for n = 1:N
        if(size(Y_p{n},1)~=0)
            max_iter = 200;
            A{n}   = N_CTF_AO_KL_sub(Y_p{n},I,F,A,lambda,opts.marg,rows{n},cols{n},n,max_iter,tol);
        end
    end
    
    %%% Update prior weights
    max_iter = 50;
    %[lambda,out_lambda] = N_CTF_AO_KL_sub_lambda(Y,F,A,lambda,opts.marg,max_iter,tol);
    [lambda,out_lambda] = N_CTF_AO_KL_sub_lambda(Y_n,F,A,lambda,marg_p,max_iter,tol);
    Out.hist_cost(iter) = out_lambda.cost;
    Out.time_instants(iter) = toc(n_tic);
    if iter>1
        if (iter == opts.max_iter ||  abs(Out.hist_cost(iter) - Out.hist_cost(iter-1))/abs(Out.hist_cost(iter-1)) < opts.tol_impr )
            Out.iter = iter;
            Out.hist_cost(iter+1:end) = [];
            Out.time_instants(iter+1:end) = [];
            break;
        end
        if mod(iter,5) == 0, fprintf('Iteration : %d cost : %d  \n', iter, Out.hist_cost(iter)); end
    end
    iter = iter + 1;
    n_tic=tic;
end
end