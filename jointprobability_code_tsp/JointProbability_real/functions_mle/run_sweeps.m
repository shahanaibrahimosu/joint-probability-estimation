function [AA,ll,Out]= run_sweeps(X,AA,ll,opts)
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
Nit=opts.max_iter;%20;       %max number of iterations
eplog=opts.tol;%1e-6;     %threshold
Nstrk=2;       %"streak" of updates beow threshold required 

%extract dimensions
FF=length(ll);
N=length(AA);
Ivec=zeros(1,N);
for n=1:N
    Ivec(n)=size(AA{n},1);
end

% get long_x

long_x=[];
for i=1:length(X)
    long_x=[long_x;X{i}(:)];
end
n_tic = tic;
Out.C_KLD=zeros(Nit,1); %for keeping track of the progress of the overall KLD
                    %between the implies estimated PMF and the empirical
                    %PMF
Out.time_instants     = zeros(Nit,1);
%Xest=lamA(ll,AA);   %the (implied) estimated PMF tensor
%C_KLD(1)=-Xemp(:)'*log(Xest(:))-entXemp;    %initial value of the KLD
Out.C_KLD(1)=Loss_Coupled(X,AA,ll,opts);
Out.time_instants(1) = toc(n_tic);
act=true(length(long_x),1); %keep track of relevant elements in long_x
nstrk=0;        %streak counter
for nit=1:Nit   %outer iterations ("sweeps")
    
    
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
    %Xest=lamA(ll,AA);   %the (implied) estimated PMF tensor
    %C_KLD(nit+1)=-Xemp(:)'*log(Xest(:))-entXemp;    %new KLD cost
    Out.time_instants(nit+1) = toc(n_tic);
    Out.C_KLD(nit+1) = Loss_Coupled(X,AA,ll,opts);
    if abs(Out.C_KLD(nit)-Out.C_KLD(nit+1))/abs(Out.C_KLD(nit))<eplog %change in log(C_KLD)<epsilon?
        nstrk=nstrk+1;  %add to streak
        if nstrk>=Nstrk, break, end
    else
        nstrk=0;
    end
    disp(['sweep ' num2str(nit)])
    n_tic=tic;
end
disp(['Completed in ' num2str(nit) ' sweeps'])
Out.C_KLD=Out.C_KLD(1:nit+1);   %truncate to actual length
end
function err = Loss_Coupled(X,A,lambda,opts)
err = 0;
A1 = {};
for i = 1 : length(X)
    A1{1}  = A{opts.marg{i}(1)}*diag(lambda);
    Xest = cpdgen([A1; A(opts.marg{i}(2:end))]);
    Xemp= X{i}(X{i}~=0);Xemp=Xemp./sum(Xemp,'all');
    Xest=Xest(X{i}~=0);Xest=Xest./sum(Xest,'all');
    kl_er = Xemp(:)'*log(Xemp(:)./Xest(:));
    err = err + kl_er;
end
err = err/length(X);
end