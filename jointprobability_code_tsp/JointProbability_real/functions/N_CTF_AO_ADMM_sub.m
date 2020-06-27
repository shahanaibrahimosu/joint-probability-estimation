function [A,U] = N_CTF_AO_ADMM_sub(G,V,A,U,n,opts,rho,max_iter,Z,W)
[ ~, k ] = size(A);
% ovewrite
% rho = trace(G)/k;
L        = chol(G + rho*eye(k), 'lower');
tol      = 1e-7;
tol1     =1e-3;
epsilon=1/max(eig(G + rho*eye(k)));
max_inner_iter=10;
for itr = 1:max_iter
    A0 = A;
    Atm = L'\ ( L\ ( V + rho*(A+U)') );
    At=Atm;%randn(size(A));
    for itr_inner = 1:max_inner_iter
        %At = At-epsilon*((W)*(W')*At-W*Z'-rho*(A+U-At));
        At = At-epsilon*(G*At-V-rho*(A+U-At));
        if(cost(Z,W,At,U,A,rho)<tol1)
            break;
        end
    end
    At;
    A  = proxr(At'-U, opts, n);
    U  = U + A - At';
    r  = A - At';
    s  = (A - A0);
    if  norm(r(:)) < tol  && norm(s(:)) < tol
        break
    end
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

function c = cost(Z,W,At,U,A,rho)
    c=norm(Z-At'*W)+rho*norm(A+U-At');
end