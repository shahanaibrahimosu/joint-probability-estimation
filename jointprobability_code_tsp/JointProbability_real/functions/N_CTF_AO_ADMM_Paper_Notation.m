function [A,lambda,Out] = N_CTF_AO_ADMM_Paper_Notation(Y,I,F,opts)
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
if ~isfield(opts,'tol_impr');   opts.tol_impr   = 1e-5;          end

rho    = opts.rho;
A      = opts.A0;          % Initial tensor factors
lambda = opts.l0;     % Prior of hidden variable

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

% cost and relative cost
Out.hist_cost         = zeros(opts.max_iter,1);
Out.hist_rel_cost     = zeros(opts.max_iter,1);
Out.time_instants     = zeros(opts.max_iter,1);

iter = 1;
rows = cell(N,1);
cols = cell(N,1);

for n = 1:N
    for i=1: size(opts.marg,1)
        [~,c]= find(opts.marg{i} == n);
        if(~isempty(c))
            rows{n}= [rows{n} i];
            cols{n}= [cols{n} c];
        end
    end
end

n_tic = tic;
while(1)
    %%%  Solve each subproblem with ADMM
    
    for n = 1:N
        G = computeG(opts.marg,rows{n},cols{n},lambda,GG);
        V = computeV(opts.marg,rows{n},cols{n},lambda,A,Y(rows{n}),I(n));
        max_iter = 500;
        [A{n}, U{n}] = N_CTF_AO_ADMM_sub_Paper_Notation(G,V,A{n},U{n},n,opts,rho(1),max_iter);
        GG{n} = A{n}'*A{n};
    end
    
    G = computeG_lambda(opts.marg,GG,F);
    V = computeV_lambda(opts.marg,A,Y,F);
    max_iter = 100;
    
    [lambda, U{N+1}] = N_CTF_AO_ADMM_sub_Paper_Notation(G,V,lambda',U{N+1},N+1,opts,rho(2),max_iter);
    lambda = lambda';

    [Out.hist_cost(iter),Out.hist_rel_cost(iter)] = Loss_Coupled(Y,A,opts,lambda,nrm_sqr);
    Out.time_instants(iter) = toc(n_tic);
    if iter>1
        if (iter == opts.max_iter ||  abs(Out.hist_rel_cost(iter) - Out.hist_rel_cost(iter-1))/abs(Out.hist_cost(iter-1)) < opts.tol_impr )
            Out.iter = iter;
            Out.time_instants(iter+1:end) = [];
            Out.hist_rel_cost(iter+1:end) = [];
            Out.hist_cost(iter+1:end)     = [];
            break;
        end
        if mod(iter,5) == 0, fprintf('Iteration : %d rel cost : %d rel cost diff : %d \n', iter, Out.hist_rel_cost(iter), abs(Out.hist_rel_cost(iter) - Out.hist_rel_cost(iter-1))); end
    end
    iter = iter + 1;
     n_tic=tic;
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
