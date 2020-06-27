function [A,A_prev] = AGD_nesterov_sub_lambda(A,n,opts,max_iter,Z,W,A_prev)
[ ~, k ] = size(A);
% ovewrite
% rho = trace(G)/k;
omega=0.3;
sigma=0.5;
beta=0.5;
tol      = 1e-7;
epsilon=1/(max(eig((W')*W)));
A_old=A;
cost_prev= cost(Z,W,A');
for itr = 1:max_iter
    A = A_old+omega*(A_old-A_prev);
    grad_A = A*(W')*(W)-Z'*W;
    A = A-epsilon*(grad_A);
    A  = proxr(A', opts, n);
    A=A';
    if(cost(Z,W,A')<tol)
        break;
    end
    cost_new = cost(Z,W,A');
    d=A_old-A;
    A_prev = A_old;
%     while((cost_prev-cost_new) < sigma* grad_A(:)*d(:)')
%         epsilon=epsilon*beta;
%         A = A-epsilon*grad_A;
%         A  = proxr(A', opts, n);
%         A=A';
%         cost_new=cost(Z,W,A');
%         d = A_old-A;
%         if(epsilon< 1e-5)
%             A = A_old;
%             break;
%         end
%     end
end
end

function A = proxr(Ab,opts,n)
switch opts.constraint{n}
    case 'nonnegative'
        A  = max(0, Ab);
    case 'simplex'
        A = reshape(ProjectOntoSimplex(Ab(:),1),size(Ab));
    case 'simplex_col'
        %         A = ProjectOntoSimplex(Ab-1e-5,1-size(Ab,1)*1e-5);
        %         A = A+1e-5;
        A = ProjectOntoSimplex(Ab,1);
end
end

function c = cost(Z,W,At)
    c=norm(Z-W*At);
end
