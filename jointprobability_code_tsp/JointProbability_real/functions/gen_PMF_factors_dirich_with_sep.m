function [A,prior] = gen_PMF_factors_dirich_with_sep(I,F,eps)
N     = length(I);
A     = cell(N,1);

prior = ones(F,1);

for n = 1 : N
    alpha = ones(1,F);
    %A{n} = transpose(drchrnd(alpha,F));
    A{n} = (drchrnd(alpha,I(n)));
    A{n}=A{n}*diag(1./sum(A{n},1));
end

prior = prior/sum(prior);
if(eps==0.1)
    epsilon=0.0005;
elseif(eps==0.2)
    epsilon=0.015;
elseif(eps==0.3)
     epsilon=0.030;
else
    epsilon=0.100;
end
%0.0005 for epsilon=0.1, 0.015 for epsilon=0.2, 0.030 for epsilon=0.3
mu = 0.4;
u=randsample([4 5],1);
A{u}(F+1:end,:) = A{u}(F+1:end,:)*diag(1./sum(A{u}(F+1:end,:),1)).*(1-(mu+4*epsilon));
L = mu*eye(F);
L(L==0)=epsilon;
A{u}(1:F,:)=L;
A{u}=A{u}*diag(1./sum(A{u},1));
%norm([1 0 0 0 0]-A{u}(1,:)/sum(A{u}(1,:)))
end
