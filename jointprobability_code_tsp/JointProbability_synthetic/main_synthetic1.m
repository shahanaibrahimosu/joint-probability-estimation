clc
clearvars
close all

addpath functions
addpath functions/tensorlab

%Test Variables
N = 15; %Number of Variables
I_ind_list = [10] ; %Number of values each variable can take
F_true = 10; %Rank of the matrix factorization
M = 3; %Parameter used for SPA/XRAY initlaization

% Number of simulations
n_sim = 20;

%No of samples
n_items_list=[1000,10000,100000, 1000000];

%Initialization Method
init_method = 1; %0- Random Initialization, % 1 - SPA initialization

%Noise addition
noise_flag = 1; %0- No noise in the marginals, 1-Noise in Marginals

%GD parameters
max_iter=1000;
tol=10^-6;

%Observed parameter
p_obs = 0.5*ones(N,1);

ii=1;
for n_items= n_items_list
for I_ind=I_ind_list
I = I_ind*ones(1,N); % Size of the Joint PMF Tensor
for s = 1:n_sim
    
    %True Parameters
    [A_true,l_true]  = gen_PMF_factors_dirich(I,F_true);    % True factors
    l_true = ones(F_true,1)/F_true;
    
    % Get the combinations of pairs
    marg    = num2cell(combnk(1:length(I),2),2);
    marg_t    = num2cell(combnk(1:length(I),3),2);
   
    %Generating marginals
    if(noise_flag==0)
        Y = get_true_marg_fast(A_true,marg,I,l_true);% Get the lower-dimensional tensors
        Y_t = get_true_marg_fast_tensor(A_true,marg_t,I,l_true);% Get the lower-dimensional tensors
    else
        %Generate data samples
        dataset = generate_dataset(A_true,l_true,n_items,I,F_true,p_obs);
        [Y,~,~] = get_obs_marg_fast(dataset,marg,I);  
        Y_t = get_obs_marg_fast(dataset,marg_t,I); 
    end
    
   %Run the initialization using singleSPA algorithm 
   opts={};
   opts.marg=marg;
   opts.max_iter=max_iter;
   opts.tol=tol;
   disp('Running M-SPA Algorithm....');
   a=tic;
   [A_est,l_est] = SingleSPA(Y,I,M,F_true,opts);
   time_mspa(s,ii)=toc(a);
    for i=1:N
       A_est{i}(isnan(A_est{i}))=0;
    end
   l_est(isnan(l_est))=0;
   %Resolve permutation ambiguity
   [rec_error1(s,ii),A_est_p,l_est_p] = getMSE_nomre(A_est,l_est,A_true,l_true);
   disp(['Mean Estimation Error - M-SPA: ',num2str(rec_error1(s,ii))]);       
  
   
     %Run the KL algorithm
    [A0,l0]= gen_PMF_factors(I,F_true);
    l0=1/F_true*ones(F_true,1);
    % Algorithm options
    opts2 = {};
    for n = 1:N
        opts2.constraint{n} = 'simplex_col';
    end
    opts2.constraint{N+1}   = 'simplex';
    opts2.max_iter = 200;
    opts2.rho = [10 10]; 
    opts2.A0=A_est;opts2.l0=l_true;
    for i=1:length(opts2.A0)
       G =opts2.A0{i};
       G = max( G, 10^-6 );
       t=sum(G,1);
       G = G*diag(1./t);
       opts2.A0{i}=G;
    end
    opts2.marg    = marg;
    opts2.computeCost=1;
    opts2.computeMSE=0;
    opts2.computeCostInterval=20;
    opts2.tol_impr = 1e-6;
    disp('Running KL_Matrix Algorithm....');
    [A_est_LS_m,l_est_LS_m,Out_LS] = N_CTF_AO_KL(Y,I,F_true,opts2);
    temp=cumsum(Out_LS.time_instants);
    time_m(s,ii)=temp(end);
   %Resolve permutation ambiguity
   [rec_error2(s,ii),A_est_LS_p,l_est_LS_p] = getMSE_nomre(A_est_LS_m,l_est_LS_m,A_true,l_true);
   disp(['Mean Estimation Error - KLMatrix: ',num2str(rec_error2(s,ii))]);

   
    %Run EM algorithm
    [A] = convert_for_comp(dataset');
    Z = cell(N,1);
    for n=1:N
        Z{n}=zeros(n_items,I(n));
    end
    for i = 1:size(A,1)
        Z{A(i,2)}(A(i,1),A(i,3)) = 1;
    end 
    opts={};
    opts.Nround = 30; %number of EM iterations
    opts.A_true=A_true;
    opts.l_true=l_true;
    opts.tol = 10^-6;
    for i=1:length(A_est)
       G =A_est{i};
       G = max( G, eps );
       t=sum(G,1);
       G = G*diag(1./t);
       A_est{i}=G;
    end
    [A_est,l_est,timestamps] = run_EM_mod(A_est,l_est,Z,opts);
    time_em(s,ii)=sum(timestamps);
    %Resolve permutation ambiguity
    [rec_error3(s,ii),A_est_p,l_est_p] = getMSE_nomre(A_est,l_est,A_true,l_true);
    disp(['Mean Estimation Error - SPA-EM: ',num2str(rec_error3(s,ii))]);
    
    
    %Run EM algorithm
    [A] = convert_for_comp(dataset');
    Z = cell(N,1);
    for n=1:N
        Z{n}=zeros(n_items,I(n));
    end
    for i = 1:size(A,1)
        Z{A(i,2)}(A(i,1),A(i,3)) = 1;
    end
    opts={};
    opts.Nround = 30; %number of EM iterations
    opts.A_true=A_true;
    opts.l_true=l_true;
    opts.tol = 10^-6;
    [A0,l0]= gen_PMF_factors(I,F_true);
    %l0=1/F_true*ones(F_true,1);
    [A_est,l_est,timestamps] = run_EM_mod(A0,l0,Z,opts);
    time_rand_em(s,ii)=sum(timestamps);
    %Resolve permutation ambiguity
    [rec_error4(s,ii),A_est_p,l_est_p] = getMSE_nomre(A_est,l_est,A_true,l_true);
    disp(['Mean Estimation Error - RAND-EM: ',num2str(rec_error4(s,ii))]);


    %Run the KLTensor algorithm
    [A0,l0]= gen_PMF_factors(I,F_true);
    l0=1/F_true*ones(F_true,1);
    % Algorithm options
    opts2 = {};
    for n = 1:N
        opts2.constraint{n} = 'simplex_col';
    end
    opts2.constraint{N+1}   = 'simplex';
    opts2.max_iter = 200;
    opts2.rho = [10 10]; 
    opts2.A0=A0;opts2.l0=l0;
    opts2.marg    = marg_t;
    opts2.computeCost=1;
    opts2.computeMSE=0;
    opts2.computeCostInterval=20;
    opts2.tol_impr = 1e-6;
   disp('Running KL_Tensor Algorithm....');
    [A_est_LS,l_est_LS,Out_LS] = N_CTF_AO_ADMM_Paper_Notation(Y_t,I,F_true,opts2);
    temp=cumsum(Out_LS.time_instants);
    time_t(s,ii)=temp(end);
    %Resolve permutation ambiguity
   [rec_error5(s,ii),A_est_LS_p,l_est_LS_p] = getMSE_nomre(A_est_LS,l_est_LS,A_true,l_true);
   disp(['Mean Estimation Error - LSTensor: ',num2str(rec_error5(s,ii))]);
%   
    opts={};
    opts.Nround = 30; %number of EM iterations
    opts.A_true=A_true;
    opts.l_true=l_true;
    opts.tol = 10^-6;
    [A_est,l_est,timestamps] = run_EM_mod(A_est_LS,l_est_LS,Z,opts);
    time_em_ctd(s,ii)=sum(timestamps);
    %Resolve permutation ambiguity
    [rec_error6(s,ii),mre_error6(s,ii),A_est_p,l_est_p] = getMSE(A_est,l_est,A_true,l_true,1);
    disp(['Mean Estimation Error - CTD-EM: ',num2str(rec_error6(s,ii))]);
    disp(['Mean Relative Error - CTD-EM: ',num2str(mre_error6(s,ii))]);
 
end
ii=ii+1;
end
end


 
    figure(3)
   loglog(n_items_list,nanmean(rec_error1,1),'-or','LineWidth',2)
    hold on;
    loglog(n_items_list,nanmean(rec_error2,1),'-*b','LineWidth',2)
    hold on;
    loglog(n_items_list,nanmean(rec_error3,1),'-dc','LineWidth',2)
    hold on;
    loglog(n_items_list,nanmean(rec_error4,1),'-<m','LineWidth',2)
    hold on;
    loglog(n_items_list,nanmean(rec_error5,1),'->g','LineWidth',2)


    grid on;
        title(['N= ',num2str(N),', I=',num2str(I(1)),', F=',num2str(F_true),', M=',num2str(M)]);
    xlabel("Number of items");
    ylabel("MSE");
    legend("CNMF-SPA","CNMF-OPT","SPA-EM","RAND-EM","CTD");
    set(gca,'fontsize',14);
