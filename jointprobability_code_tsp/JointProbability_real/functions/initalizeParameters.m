function [A0,l0] = initializeParametrs(I,F,Y,marg)

    X_cell =get_golden_matrix(Y,marg,M,I);
    X=X_cell{1,1};
    X_sum=sum(X,1);
    
       if(XRAY~=1)
        X = X*diag(1./sum(X,1));
        [l,jk,kl]=FastSepNMF(X,F_true,0);
        W_est= X(:,l);
        H_est=(pinv(W_est)*X);
        H_est = H_est';
        H_est = diag(transpose(X_sum))*H_est;

    else        
        [l,H_est]=FastConicalHull(X,F_true);
        W_est= X(:,l); 
        H_est = H_est';        
    end
    

        
    A_est_G = cell(N,1);
    for n=1:M
        A_est_G{n} = W_est(sum(I(1:n-1))+1:sum(I(1:n)),:);
        A_est_G{n} = A_est_G{n}*diag(1./sum(A_est_G{n},1));
    end

    for n=M+1:N
        A_est_G{n} = H_est(sum(I(M+1:n-1))+1:sum(I(M+1:n)),:);
        A_est_G{n} = A_est_G{n}*diag(1./sum(A_est_G{n},1));
    end