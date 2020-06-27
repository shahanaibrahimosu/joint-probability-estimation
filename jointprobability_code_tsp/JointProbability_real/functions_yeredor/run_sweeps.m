%run_sweeps
%Alternating directions miniization with respect to the factor matrices and
%loading vector.
%
%For minimizing wrt An we need to minimize a function of the form 
%            -sum_{f=1}^FF long_x'*(long_HF(:,:,f)*An(:,f)
%subject to the simplex constraint on An, where long_HF(Nn,In,FF) is 
%constructed from the other (fixed) matrices and the loading vector. 
%Here Nn is the number of elements in long_x which are relevant to the 
%minimization wrt An %(namely elements coming from (j,k,l) triplets in 
%which either j, k or l equal n. The dimensions of An are In by FF.
%
%For minimizing wrt ll we need to minimize a function of the form 
%            -long_x'*(long_HF*ll)
%subject to the simplex constraint on ll, where long_HF(NN,FF) is 
%constructed from the (fixed) factor matrices. 
%Here NN is the entire length of long_x, since all elements are relevant
%to the minimization wrt ll.
%
%This code sweeps over all factor matrices and the loading vector (each at 
%a time, with the others fixed) and constructs the HF matrices for the
%inner minimizations (done by fmyminKLD)
%
%By Arie Yeredor, Nov. 2018, revised May 2020

%parameters for the outer minimization
Nit=1000;       %max number of iterations
eplog=1e-5;     %threshold
Nstrk=10;       %"streak" of updates beow threshold required 

C_KLD=zeros(Nit,1); %for keeping track of the progress of the overall KLD
                    %between the implies estimated PMF and the empirical
                    %PMF
Xest=lamA(ll,AA);   %the (implied) estimated PMF tensor
C_KLD(1)=-Xemp(:)'*log(Xest(:))-entXemp;    %initial value of the KLD

act=true(length(long_x),1); %keep track of relevant elements in long_x
nstrk=0;        %streak counter
for nit=1:Nit   %outer iterations ("sweeps")
    disp(['sweep ' num2str(nit)])
    
    %loop over each of the N factor matrices {An}:
    for nmat=1:N 
        In=Ivec(nmat);  %the dimension of An
        I=eye(In);
        long_HF=[];
        actn=act;   %to define "active" indices of long_x, to which this matrix is relevant
        base=0;
        for j=1:N-2
            for k=j+1:N-1
                for l=k+1:N
                    Ijkl=prod(Ivec([j k l])); %number of elements of the respective triplet PMF
                    HF=zeros(Ijkl,In,FF);
                    switch nmat
                        case j
                            for f=1:FF
                                HF(:,:,f)=ll(f)*kron(AA{l}(:,f),kron(AA{k}(:,f),I));
                            end
                            long_HF=vertcat(long_HF,HF);
                        case k
                            for f=1:FF
                                HF(:,:,f)=ll(f)*kron(AA{l}(:,f),kron(I,AA{j}(:,f)));
                            end
                            long_HF=vertcat(long_HF,HF);
                        case l
                            for f=1:FF
                                HF(:,:,f)=ll(f)*kron(I,kron(AA{k}(:,f),AA{j}(:,f)));
                            end
                            long_HF=vertcat(long_HF,HF);
                        otherwise
                            %this triplet is irrelevant to An, mark the respective 
                            %elements of long_x as irrelevant
                            actn(base+[1:Ijkl])=false; %not "active"
                    end
                    base=base+Ijkl;
                end
            end
        end
        AA{nmat}=fmyminKLD(long_x(actn),long_HF);   %find the minimizing An
    end
    
    %now minimize wrt ll
    %no need for actn here, since all elements (full long_x) participate
    long_HF=[];
    for j=1:N-2
        for k=j+1:N-1
            for l=k+1:N
                Ijkl=prod(Ivec([j k l]));
                HF=zeros(Ijkl,FF);
                for f=1:FF
                    HF(:,f)=kron(AA{l}(:,f),kron(AA{k}(:,f),AA{j}(:,f)));
                end
                long_HF=vertcat(long_HF,HF);
            end
        end
    end
    ll=fmyminKLD(long_x,long_HF);   %find the minimizing ll
    
    %check for the stopping criterion, monitoring C_KLD
    Xest=lamA(ll,AA);   %the (implied) estimated PMF tensor
    C_KLD(nit+1)=-Xemp(:)'*log(Xest(:))-entXemp;    %new KLD cost
    if log10(C_KLD(nit))-log10(C_KLD(nit+1))<eplog %change in log(C_KLD)<epsilon?
        nstrk=nstrk+1;  %add to streak
        if nstrk>=Nstrk, break, end
    else
        nstrk=0;
    end
end
disp(['Completed in ' num2str(nit) ' sweeps'])
C_KLD=C_KLD(1:nit+1);   %truncate to actual length