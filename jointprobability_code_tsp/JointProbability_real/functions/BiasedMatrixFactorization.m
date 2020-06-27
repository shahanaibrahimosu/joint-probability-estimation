function [A,B,cost_it,RMSE3,RMSE2,MAE2] = BiasedMatrixFactorization(G,r1,c1,P3,r3,c3,P4,r4,c4,lbd,F)
%COLAB_FILTERING implements a colaborative filtering algorithm with missing
%values, the algorithm solves min ||X-AB'|| + lbd*(||A||+||B||), the least
%squares are calculated only the the position where observation is
%available. Solving the problem with Alternating Least Squares(ALS).
%   X: data matrix, with missing values represented by NaN
%   F: the length of the latent factors.

% r1,c1 are the positions of missing values
% lbd is the regularization parameter
% G is the input matrix

X = G;

avG = sum(X(X~=0))/nnz(X);
for i=1:length(r1);
    X(r1(i),c1(i)) = NaN;
end

mask = isnan(X); % missing entries of X

X0 = X;
X0(mask) = avG;       % fill missing entries with average
mask = ~mask;         % now mask have the nz entries

tot = sum(sum(mask)); % total number of observed entries
[m,n] = size(X);

A = zeros(m,F);       % initialize factors
B = randn(n,F);

it_tot = 100;
it = 1;

% User and item bias
a = zeros(1,m)';
b = zeros(1,n)';
u = 0;
test = inf;
cost_it = zeros(1,it_tot);

while it < it_tot
    it;
    it  = it+1;
    %% Solving for A: user
    for i = 1:m
        Bi = B(mask(i,:),:);
        Xi = X(i,mask(i,:))'-u-a(i)-b(mask(i,:));
        XX = [Bi; lbd*eye(F)];
        yy = [Xi; zeros(F,1)];
        A(i,:) = (XX\yy)';
    end
    %% Solving for B: item
    for i = 1:n
        Ai = A(mask(:,i),:);
        Xi = X(mask(:,i),i)-u-b(i)-a(mask(:,i));
        XX = [Ai; lbd*eye(F)];
        yy = [Xi; zeros(F,1)];
        B(i,:) = (XX\yy)';
    end
    %% Update overall average
    u = sum(sum(mask.*(X0-A*B'-a*ones(1,n)-ones(m,1)*b')))/tot;
    
    %% Update user bias
    for i=1:m
        num = sum(mask(i,:));
        Xi = X(i,mask(i,:))'-u-b(mask(i,:))-B(mask(i,:),:)*A(i,:)';
        yy = [Xi;zeros(num,1)];
        XX = [ones(num,1);lbd*ones(num,1)];
        a(i) = XX\yy;
    end
    %% Update item bias
    for i = 1:n
        num = sum(mask(:,i));
        Xi = X(mask(:,i),i)-u-a(mask(:,i))-A(mask(:,i),:)*B(i,:)';
        yy = [Xi;zeros(num,1)];
        XX = [ones(num,1);lbd*ones(num,1)];
        b(i) = XX\yy;
    end
    
    %% Calculating residuals
    cost = norm(mask.*(X0-u-a*ones(1,n)-ones(m,1)*b'-A*B'),'fro')^2 + ...
        lbd^2*(norm(A,'fro')^2+norm(B,'fro')^2+norm(a)^2+norm(b)^2);
    cost_it(it) = cost;
    
    %% calculate RMSE
    Bt = B';
    P3p = [];
    for i = 1:length(r3);
        P3p = [P3p A(r3(i),:)*Bt(:,c3(i))+a(r3(i))+b(c3(i))+u];
    end
    
    Err = 0;
    count = 0;
    for i = 1:length(P3)
        if any(G(r3(i),:)~=0)
            Err = Err + (P3p(i) - P3(i))^2;
            count = count + 1;
        end
    end
    RMSE = sqrt(Err/count);
    test = [test, RMSE];
    
   
    if it>1
        if (test(it-1) - test(it)) < 0
            RMSE3 = RMSE;
            break
        end
    end
end


%% P4
P4p = [];
for i = 1:length(r4);
    P4p = [P4p; A(r4(i),:)*Bt(:,c4(i))+a(r4(i))+b(c4(i))+u];
end

% ratings between 1 and 5
P4p(P4p>5)   = 5;
P4p(P4p<1)   = 1;

%% RMSE
count = 0;
Err = 0;
for i = 1:length(P4);
    if any(G(r4(i),:)~=0)
        Err = Err + (P4p(i) - P4(i))^2;
        count = count + 1;
    end
end
RMSE2 = sqrt(Err/count)

%% MAE
count = 0;
Err2 = 0;
for i = 1:length(P4);
    if any(G(r4(i),:)~=0)
        Err2 = Err2 + abs(P4p(i) - P4(i));
        count = count + 1;
    end
end
MAE2 = Err2/count

count;
end