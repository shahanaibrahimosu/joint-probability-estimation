function V = computeV_lambda(marg,A,Y,F)
V = zeros(F,1);
for i = 1:size(marg,1)
    V = V + khatrirao(A(marg{i}),'r')'*Y{i}(:);
end

% V = zeros(F,1);
% for i = 1:size(marg,1)
%     for f=1:F
%         vect = cell(size(marg,2),1);
%         for j = 1:size(marg,2)
%             vect{j} = A{marg(i,j)}(:,f)';
%         end
%         V(f) = V(f) + tmprod(Y{i},vect,1:size(marg,2));
%     end
% end