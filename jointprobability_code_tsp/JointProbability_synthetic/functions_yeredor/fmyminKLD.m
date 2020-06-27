function AA=fmyminKLD(xx,HH);
%inputs: xx HH
%output: AA
%
%xx is an NNx1 vector
%HH is an NNxIIxFF array (of FF matrices);
%this function finds an IIxFF matrix AA (vector if FF=1) minimizing
%-sum_{f=1}^FF xx'*(HH(:,:,f)*AA(:,f)
%where AA is subject to the probability simplex constraint, namely all of
%its elements are non-negative, 1^T*AA=1.
%
%See the paper 
%A. Yeredor and M. Haardt, "Estimation of a Low-Rank Probability-Tensor 
%from Sample Sub-Tensors via Joint Factorization Minimizing the 
%Kullback-Leibler Divergence", EUSIPCO 2019
%for details of the minimization process
%
%By Arie Yeredor, Nov. 2018, revised May 2020

%some parameters:
mNit=1000;  %maximum number of iterations
alpha=0.01; %gradient step size factor
myep=1e-7;  %update threshold for stopping

[NN II FF]=size(HH);
hbar=sum(HH(:,1,:),3);
Q=toeplitz([1 -1 zeros(1,II-2)],[1 zeros(1,II-2)]);
Hbar=zeros(NN,II-1,FF);
for f=1:FF
    Hbar(:,:,f)=HH(:,:,f)*Q;
end

BB=repmat([II:-1:0]'/II,1,FF);
B=BB(2:II,:);
oldB=B;
G=zeros(II-1,FF);
for mnit=1:mNit
    HB=zeros(NN,1);
    for f=1:FF
        HB=HB+Hbar(:,:,f)*B(:,f);
    end
    for f=1:FF
        G(:,f)=(xx./(hbar-HB))'*Hbar(:,:,f);
    end
    [Gs,Is]=sort(G*alpha,1,'descend');
    for f=1:FF
        bb=BB(:,f);
        for ii=1:II-1
            stp=Gs(ii,f);
            bix=Is(ii,f);
            bb(bix+1)=bb(bix+1)-stp;
            bb(bix+1)=max(min(bb(bix+1),bb(bix)-myep),bb(bix+2)+myep);
        end
        BB(:,f)=bb;
    end
    B=BB(2:II,:);
    DB=B-oldB;
    if max(abs(DB))<myep
        break
    end
    oldB=B;
end
AA=[1;zeros(II-1,1)]-Q*B;
    
    
    



