function [A,U] = N_CTF_AO_ADMM_lambda(G,V,A,U,n,opts,rho,max_iter,Z,W)
[ ~, k ] = size(A);
% ovewrite
% rho = trace(G)/k;
L        = chol((W')*W + rho*eye(k), 'lower');
tol      = 1e-7;
for itr = 1:max_iter
    A0 = A;
    At = L'\ ( L\ ( Z'*W + rho*(A+U)') );
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
