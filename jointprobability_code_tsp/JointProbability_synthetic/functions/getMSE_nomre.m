function [rec_error,A,l] = getMSE(A,l,A_true,l_true)
N=size(A,1);
F=size(A{1},2);

%Resolve permutation ambiguity
A_est_cat   = concatenate(A);
A_true_cat  = concatenate(A_true);
[~,Pm] = perm2match(A_est_cat,A_true_cat);
for n=1:N
    A{n} = A{n}*Pm;
end
l = Pm*l;

 A_est=A;
% for i=1:length(A_est)
%    G =A_est{i};
%    G = max( G, 10^-6);
%    t=sum(G,1);
%    G = G*diag(1./t);
%    A_est{i}=G;
% end

% %% compute error 1
% A_true{1}=A_true{1}*diag(l_true);
% T_true=cpdgen(A_true);
% A_est{1}=A_est{1}*diag(l);
% T_est=cpdgen(A_est);
% rec_error1 = frob(T_true-T_est)/frob(T_est);
%rec_error1 = sum(tens2vec(T_true).*log(tens2vec(T_true)./tens2vec(T_est)));




rec_error=norm(l-l_true)^2;
% Compute error
for n=1:N
    rel_factor_error = (norm(A{n}(:) - A_true{n}(:))^2)/F;
    rec_error =  rec_error + rel_factor_error/N;
end

% % Compute error1
% A_true{1}=A_true{1}*diag(l_true);
% T_true=cpdgen(A_true);
% A{1}=A{1}*diag(l);
% T_est=cpdgen(A);
% rec_error1 = sum(tens2vec(T_true).*log(tens2vec(T_true)./tens2vec(T_est)));
% 
% rec_error2=sum(l_true.*log(l_true./l));
% % Compute error
% for n=1:N
%     for f=1:F
%         rel_factor_error = sum(A_true{n}(:,f).*log(A_true{n}(:,f)./A{n}(:,f)));
%         rec_error2 =  rec_error2 + rel_factor_error;
%     end
% end
% rec_error2=rec_error2/(N*F);



end