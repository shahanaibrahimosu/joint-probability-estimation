function [A,prior] = gen_PMF_factors_dirich(I,F)
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
% epsilon=0;
% mu = 0.4;
% u=randsample([4 5],1);
% A{u}(F+1:end,:) = A{u}(F+1:end,:)*diag(1./sum(A{u}(F+1:end,:),1)).*(1-(mu+4*epsilon));
% L = mu*eye(F);
% L(L==0)=epsilon;
% A{u}(1:F,:)=L;
% A{u}=A{u}*diag(1./sum(A{u},1));
end
