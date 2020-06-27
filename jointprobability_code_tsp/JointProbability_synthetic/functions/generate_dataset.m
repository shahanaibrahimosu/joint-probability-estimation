function [dataset,dataset_ML]=generate_dataset(A_true,l_true,n_items,I,F_true,p_obs)
    N=length(I);
    dataset = zeros(n_items,N+1);
    dataset(:,N+1) = sum(repmat(rand(1,n_items),F_true,1) > repmat(cumsum(l_true),1,n_items),1)+1;  
    f_data = dataset(:,N+1);
    for n=1:N
        for i=1:F_true
            ind = dataset(:,N+1)==i;
            s_ind = sum(ind);
            dataset(ind,n) = sum(repmat(rand(1,s_ind),I(n),1)> repmat(cumsum(A_true{n}(:,i)),1,s_ind),1)+1;
        end
    end 
    dataset = dataset(:,1:N);% Hide the last column (hidden variable)
    f=dataset';
    for i=1:N
        mask = binornd(1,p_obs(i),1,n_items);
        f(i,:) = f(i,:).*mask;
    end
    dataset = f';
    dataset_ML = [dataset f_data];
%     f=dataset_ML';
%     for i=1:N
%         mask = binornd(1,p_obs(i),1,n_items);
%         f(i,:) = f(i,:).*mask;
%     end
%     dataset_ML = f';


end