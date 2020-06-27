function [RMSE2,MAE2] = complete_rating_bmf(A,B,a,b,u,G,P4,r4,c4)


%% P4
Bt = B';
P4p = [];
for i = 1:length(r4)
    P4p = [P4p; A(r4(i),:)*Bt(:,c4(i))+a(r4(i))+b(c4(i))+u];
end

% ratings between 1 and 5
P4p(P4p>5)   = 5;
P4p(P4p<1)   = 1;

%% RMSE
count = 0;
Err = 0;
for i = 1:length(P4)
    if any(G(r4(i),:)~=0)
        Err = Err + (P4p(i) - P4(i))^2;
        count = count + 1;
    end
end
RMSE2 = sqrt(Err/count);

%% MAE
count = 0;
Err2 = 0;
for i = 1:length(P4)
    if any(G(r4(i),:)~=0)
        Err2 = Err2 + abs(P4p(i) - P4(i));
        count = count + 1;
    end
end
MAE2 = Err2/count;

count;
end