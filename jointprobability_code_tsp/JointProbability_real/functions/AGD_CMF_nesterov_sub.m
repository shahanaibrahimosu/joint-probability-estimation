function [A_new,A_prev] = AGD_nesterov_sub(A,n,opts,max_iter,Z,W,l,XX,A_prev)
[ ~, k ] = size(A);
% ovewrite
% rho = trace(G)/k;
omega=0.5;
sigma=0.5;
beta=0.5;
D=diag(l);
tol      = 1e-7;
epsilon=1/(max(eig((W)*(W'))));
cost_prev = cost(Z,W,A,D,XX);
A_old=A;
for itr = 1:max_iter
    A = A_old+omega*(A_old-A_prev);
    grad_A=A*(W)*(W')-Z*W';%-2*XX*D*A+2*A*D*A'*D*A;
    A_new = A-epsilon*grad_A;
    A_new = proxr(A_new, opts, n);
    cost_new = cost(Z,W,A_new,D,XX);
    d=A_old-A_new;
    A_prev = A_old;
%     while((cost_prev-cost_new) < sigma* grad_A(:)'*d(:))
%         epsilon=epsilon*beta;
%         A_new = A-epsilon*grad_A;
%         A  = proxr(A, opts, n);
%         cost_new=cost(Z,W,A,D,XX);
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

function c = cost(Z,W,At,D,XX)
    c=norm(Z-At*W);%+norm(XX-At'*D*At);
end