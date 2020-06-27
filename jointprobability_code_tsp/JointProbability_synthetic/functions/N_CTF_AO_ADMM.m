function [A,lambda,Out] = N_CTF_AO_ADMM(Y,I,F,opts,A_true)
% input:
%       y   : measurement vector
%       I   : number of rows of each factor
%       F   : rank
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

%% Parameters
N  = length(I);
U  = cell(N+1,1); % scaled dual variables
GG = cell(N,1);

if ~isfield(opts,'max_iter'),   opts.max_iter   = 1000;          end
if ~isfield(opts,'tol_impr');   opts.tol_impr   = 1e-7;          end

rho    = opts.rho;
A      = opts.A0;          % Initial tensor factors
lambda = opts.l0;     % Prior of hidden variable
n_iter = opts.max_iter;
tol    = opts.tol_impr;

nrm_sqr = 0;
for i = 1 : size(opts.marg,1)
    nrm_sqr = nrm_sqr + frob(Y{i},'squared');
end

for n = 1:N
    GG{n} = A{n}'*A{n};    % A^T*A cache
end

for n = 1:N
    U{n} = zeros(size(A{n}));
end
U{N+1}   = zeros(size(lambda'));



rows = cell(N,1);
cols = cell(N,1);

for n = 1:N
    for i=1: size(opts.marg,1)
        [~,c]= find(opts.marg{i} == n);
        if(~isempty(c)&& rank(Y{i})> 1)
            rows{n}= [rows{n} i];
            cols{n}= [cols{n} c];
        end
    end
end

% cost and relative cost
cost_iter             = floor(n_iter/opts.computeCostInterval);
Out.hist_cost         = zeros(cost_iter+1,1);
Out.rel_hist_cost         = zeros(cost_iter+1,1);
Out.time_stamps       = zeros(cost_iter+1,1);
Out.MSE               = zeros(cost_iter+1,1);

a = tic;
iter_c=1;
for iter=1:n_iter
    %Initial Objective Value
    if(opts.computeCost && iter==1)
        b = toc(a);
        Out.time_stamps(iter_c) = b;
        [Out.hist_cost(iter_c),Out.rel_hist_cost(iter_c)] = Loss_Coupled(Y,A,opts,lambda,nrm_sqr);
        fprintf('Iteration : %d cost : %d  \n', iter-1, Out.hist_cost(iter_c));
        if(opts.computeMSE)
            Out.MSE(iter_c) = getMSE(A,A_true); 
            fprintf('Iteration : %d MSE : %d  \n', iter, Out.MSE(iter_c));
        end
        iter_c=iter_c+1;
        a=tic;
    end  
    %%%  Solve each subproblem with ADMM

    for n = 1:N
        G = computeG(opts.marg,rows{n},cols{n},lambda,GG);
        V = computeV(opts.marg,rows{n},cols{n},lambda,A,Y(rows{n}),I(n));
        [Z,W]=computeZW(opts.marg,rows{n},cols{n},lambda,Y,A);
        max_iter = 3;
        [A{n}, U{n}] = N_CTF_AO_ADMM_sub(G,V,A{n},U{n},n,opts,rho(1),max_iter,Z,W);
        GG{n} = A{n}'*A{n};
    end
    
    G = computeG_lambda(opts.marg,GG,F);
    V = computeV_lambda(opts.marg,A,Y,F);
    [Z_l,W_l]=computeZW_lambda(opts.marg,Y,A);
    max_iter = 1;
    
    [lambda, U{N+1}] = N_CTF_AO_ADMM_sub_lambda(G,V,lambda',U{N+1},N+1,opts,rho(2),max_iter,Z_l,W_l);
    lambda = lambda';
    
   %Calculate the objective value
    if(opts.computeCost && mod(iter,opts.computeCostInterval)==0)
        b=toc(a);
        Out.time_stamps(iter_c) = b;
        [Out.hist_cost(iter_c),Out.rel_hist_cost(iter_c)] = Loss_Coupled(Y,A,opts,lambda,nrm_sqr);
        fprintf('Iteration : %d cost : %d  \n', iter, Out.hist_cost(iter_c)); 
        if(opts.computeMSE)
            Out.MSE(iter_c) = getMSE(A,A_true); 
            fprintf('Iteration : %d MSE : %d  \n', iter, Out.MSE(iter_c));
        end
        %Out.hist_cost(iter) = out_lambda.cost;
        if iter_c>1
            if (iter == n_iter ||  abs(Out.rel_hist_cost(iter_c) - Out.rel_hist_cost(iter_c-1))/abs(Out.rel_hist_cost(iter_c-1)) < tol )
                Out.iter = iter;
                Out.hist_cost(iter_c+1:end) = [];
                Out.time_stamps(iter_c+1:end)=[];
                if opts.computeMSE
                    Out.MSE(iter_c+1:end)=[]; 
                end
                break;
            end
        end        
        iter_c=iter_c+1;
        a=tic;
    end
end

end



function [err,rel_error] = Loss_Coupled(Y,A,opts,lambda,nrm_sqr)
err = 0;
A1 = {};
for i = 1 : size(opts.marg,1)
    A1{1}  = A{opts.marg{i}(1)}*diag(lambda);
    fro_er = frob(cpdres(Y{i} , [A1; A(opts.marg{i}(2:end))]));
    err = err + fro_er^2;
end
rel_error = sqrt(err/nrm_sqr);
err = 1/2*err;
end