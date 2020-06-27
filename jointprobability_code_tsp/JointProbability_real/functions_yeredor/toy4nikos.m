%toy4nikos
%By Arie Yeredor, May 2020

%setup parameters
simsetup

%calculate the true tensor and its entropy
Xtrue=lamA(lam,Aall);
findnzt=find(Xtrue(:)>0); %find non-zero elements (to avoid log(0) in entropy calculation)
entXtrue=-Xtrue(findnzt)'*log(Xtrue(findnzt));

%draw samples directly into an empirical histogram according to Xtrue
rng(1234)       %for repeatability; comment out for "true" randomness
r=rand(1,T);    %for the histogram the vectors can simply be viewed as random variables
Xc=cumsum(Xtrue(:));    %with this cumulative PMF
histin=zeros(prod(Ivec),1);
for t=1:T
    inX=find(r(t)<Xc,1);    %find its position in the histogram
    histin(inX)=histin(inX)+1;
end
elongX=histin/T;            %normalize
Xemp=reshape(elongX,Ivec);  %rehsape the histogram as a PMF tensor
%get its entropy
findnze=find(Xemp(:)>0); %find non-zero elements (to avoid log(0) in entropy calculation)
entXemp=-Xemp(findnze)'*log(Xemp(findnze));

get_triplets    %extract all triplets-PMFs, vectorize and concatenate in a long vector
FF=F;           %the presumed low rank is the true low rank (can be changed)
init_guess      %get an initial guess of the loading vector and factor matrices
run_sweeps      %rund the iterative alternating directions minimization

%measure the estimation errors (both RMS and KLD) wrt the true PMF
Xest=lamA(ll,AA);   %get the resulting estimated PMF tensor
D2true_RMS=sqrt(mean((Xtrue(:)-Xest(:)).^2));
D2true_KLD=-Xtrue(findnzt)'*log(Xest(findnzt))-entXtrue;
fprintf('Estimation error: %2.7f (RMS)   %2.7f (KLD)\n',D2true_RMS,D2true_KLD)