function [mu,d] = run_EM(mu,d,Z,Nround)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n=size(Z{1},1);
k=size(mu{1},2);
m=length(Z);
% EM update
for iter = 1:Nround
    q = zeros(n,k);
    q1 = zeros(k,k);
    for i=1:m
        tmp = mu{i};
        tmp(find(tmp ==0)) = eps;
        tmp(find(isnan(tmp))) = eps;
        tmp1 = d;
        tmp1(find(tmp1 ==0)) = eps;
        tmp1(find(isnan(tmp1))) = eps;
        q = q + Z{i}*log(tmp);
        %q = q + Z(:,:,i)*log(mu(:,:,i));
    end
    q=q+repmat(log(tmp1'),n,1);
    q = exp(q);
    q = bsxfun(@rdivide,q,sum(q,2));
    


    for i = 1:m
        mu{i} = (Z{i})'*q;
        
        mu{i} = AggregateCFG(mu{i},0);
        tmp_vec = sum(mu{i});
        indx = find(tmp_vec > 0);
        mu{i}(:,indx) = bsxfun(@rdivide,mu{i}(:,indx),tmp_vec(indx));
    end
    d = sum(q,1)';
    tmp_vec = sum(d);
    indx = find(tmp_vec > 0);
    d = bsxfun(@rdivide,d,tmp_vec);
    
%     cost(iter) = sum(q*log(d),"all");
%     for i=1:m
%         cost(iter)=cost(iter)+sum(q*log(mu{i}),"all");
%     end
end
%error_EM_predict = mean(y(valid_index) ~= (J(valid_index))');
%EM_err = error_EM_predict(end);

end

