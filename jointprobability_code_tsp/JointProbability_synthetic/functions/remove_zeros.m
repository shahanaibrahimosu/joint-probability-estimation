function [A,l]=remove_zeros(A,l)
    for i=1:length(A)
       G =A{i};
       G = max( G, eps);
       t=sum(G,1);
       G = G*diag(1./t);
       A{i}=G;
    end
    g = l;
    g = max(g,eps);
    t=sum(g);
    g = g./t;
    l=g;
end