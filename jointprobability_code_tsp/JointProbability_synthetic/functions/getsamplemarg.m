function [X1,sel1,k1,l1,J_sel1] = getsamplemarg(j,Y,I,marg,order_j,J)
    marg_array=cell2mat(marg);
    J_sel1=datasample(1:J(j));
    sel1=marg_array(order_j{j}(J_sel1),:);
    ind=sel1(sel1~=j);
    k1=ind(1);
    l1=ind(2);
    X1=Y{marg_array(:,1)==sel1(1) & marg_array(:,2)==sel1(2) & marg_array(:,3)==sel1(3)};
    
    J_sel2=datasample(1:J(j));
    sel2=marg_array(order_j{j}(J_sel2),:);
    ind=sel2(sel2~=j);
    k2=ind(1);
    l2=ind(2);
    X2=Y{marg_array(:,1)==sel2(1) & marg_array(:,2)==sel2(2) & marg_array(:,3)==sel2(3)};
endX1,sel1,k1,l1,J_sel1