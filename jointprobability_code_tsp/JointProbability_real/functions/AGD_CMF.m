function [A,lambda,Out] = AGD_CMF(Y,I,F,opts,A_true)
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


for n = 1:N
    XX{n} = A{n}*diag(lambda)*A{n}';
end

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



% cost and relative cost
cost_iter             = floor(n_iter/opts.computeCostInterval);
Out.hist_cost         = zeros(cost_iter+1,1);
Out.rel_hist_cost         = zeros(cost_iter+1,1);
Out.time_stamps       = zeros(cost_iter+1,1);
Out.MSE               = zeros(cost_iter+1,1);

a = tic;
iter_c=1;
A_prev=A;
lambda_prev=lambda;
for iter=1:n_iter
    %Initial Objective Value
    if(opts.computeCost && iter==1)
        b = toc(a);
        Out.time_stamps(iter_c) = b;
        [Out.hist_cost(iter_c),Out.rel_hist_cost(iter_c)] = Loss_Coupled(Y,A,opts,lambda,XX);
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
        [Z,W]=computeZW_CMF(opts.marg,rows{n},cols{n},lambda,Y,A);
        max_iter = 1;
        [A{n},A_prev{n}] = AGD_CMF_nesterov_sub(A{n},n,opts,max_iter,Z,W,lambda,XX{n},A_prev{n});
    end
    
    [Z_l,W_l]=computeZW_CMF_lambda(opts.marg,Y,A,XX);
    max_iter = 1;
    
     [lambda,lambda_prev] = AGD_CMF_nesterov_sub_lambda(lambda',N+1,opts,max_iter,Z_l,W_l,lambda_prev');
     lambda = lambda';
     lambda_prev=lambda_prev';
    for n=1:N
        XX{n}=A{n}*diag(lambda)*A{n}';
    end
    
   %Calculate the objective value
    if(opts.computeCost && mod(iter,opts.computeCostInterval)==0)
        b=toc(a);
        Out.time_stamps(iter_c) = b;
        [Out.hist_cost(iter_c),Out.rel_hist_cost(iter_c)] = Loss_Coupled(Y,A,opts,lambda,XX);
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



function [err,rel_error] = Loss_Coupled(Y,A,opts,lambda,XX)
err = 0;
A1 = {};
N=length(A);
for i = 1 : size(opts.marg,1)
    A1{1}  = A{opts.marg{i}(1)}*diag(lambda);
    fro_er = frob(cpdres(Y{i} , [A1; A(opts.marg{i}(2:end))]));
    err = err + fro_er^2;
end
% for i=1:N
%     err =err + norm(XX{i}-A{i}*diag(lambda)*A{i}')^2; 
% end
rel_error = sqrt(err);
err = 1/2*err;
end