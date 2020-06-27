clc
clearvars
close all

addpath functions
addpath functions/tensorlab
addpath functions_mle

%Test Variables
N = 15; %Number of Variables
I_ind_list = [10] ; %Number of values each variable can take
F_true = 10; %Rank of the matrix factorization
M = 3; %Parameter used for SPA/XRAY initlaization

% Number of simulations
n_sim = 10;

%No of samples
n_items_list=[100000, 1000000];
%n_items_list=[1000,10000];

%Initialization Method
init_method = 1; %0- Random Initialization, % 1 - SPA initialization

%Noise addition
noise_flag = 1; %0- No noise in the marginals, 1-Noise in Marginals

%GD parameters
max_iter=1000;
tol=10^-6;

%Observed parameter
p_obs = 0.50*ones(N,1);


ii=1;
for n_items= n_items_list
for I_ind=I_ind_list
I = I_ind*ones(1,N); % Size of the Joint PMF Tensor
for s = 1:n_sim
    %% Setting up
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
        [dataset,dataset_ML] = generate_dataset(A_true,l_true,n_items,I,F_true,p_obs);
        a=tic;
        [Y,~,~] = get_obs_marg_fast(dataset,marg,I);  
        b=toc(a)
        a=tic;
        Y_t = get_obs_marg_fast(dataset,marg_t,I); 
        b=toc(a)
    end
    
    %Random Init
    [A0,l0]= gen_PMF_factors(I,F_true);
   
   %% CNMF-SPA
   %Run the initialization using singleSPA algorithm 
   opts={};
   opts.marg=marg;
   opts.max_iter=max_iter;
   opts.tol=tol;
   disp('Running CNMF-SPA Algorithm....');
   a=tic;
   [A_est,l_est] = SingleSPA(Y,I,M,F_true,opts);
   time_mspa(s,ii)=toc(a);

   
   %Resolve permutation ambiguity
   [rec_error1(s,ii),A_est_p,l_est_p] = getMSE_nomre(A_est,l_est,A_true,l_true);
   disp(['Mean Estimation Error - CNMF-SPA: ',num2str(rec_error1(s,ii))]);       
   %disp(['Mean Relative Error - CNMF-SPA: ',num2str(mre_error1(s,ii))]);       
  
%    %% CNMF-KL
%      %Run the KL algorithm
%     [A0,l0]= gen_PMF_factors(I,F_true);
%     l0=1/F_true*ones(F_true,1);
%     % Algorithm options
%     opts2 = {};
%     for n = 1:N
%         opts2.constraint{n} = 'simplex_col';
%     end
%     opts2.constraint{N+1}   = 'simplex';
%     opts2.max_iter = 200;
%     opts2.rho = [10 10]; 
%     opts2.A0=A_est;opts2.l0=l_true;
%     for i=1:length(opts2.A0)
%        G =opts2.A0{i};
%        G = max( G, 10^-6);
%        t=sum(G,1);
%        G = G*diag(1./t);
%        opts2.A0{i}=G;
%     end
%     opts2.marg    = marg;
%     opts2.computeCost=0;
%     opts2.computeMSE=0;
%     opts2.computeCostInterval=20;
%     opts2.tol_impr = 1e-6;
%     disp('Running CNMF-OPT-KL Algorithm....');
%     [A_est_KL_m,l_est_KL_m,Out_KL] = N_CTF_AO_KL(Y,I,F_true,opts2);
%     temp=cumsum(Out_KL.time_instants);
%     time_m(s,ii)=temp(end);
% 
%    %Resolve permutation ambiguity
%    [rec_error2(s,ii),A_est_KL_p,l_est_KL_p] = getMSE_nomre(A_est_KL_m,l_est_KL_m,A_true,l_true);
%    disp(['Mean Estimation Error - CNMF-OPT-KL ',num2str(rec_error2(s,ii))]);
%    %disp(['Mean Relative Error - CNMF-OPT-KL: ',num2str(mre_error2(s,ii))]);
% 
%    %% CNMF-EM
%     %Run EM algorithm
%     [A] = convert_for_comp(dataset');
%     Z = cell(N,1);
%     for n=1:N
%         Z{n}=zeros(n_items,I(n));
%     end
%     for i = 1:size(A,1)
%         Z{A(i,2)}(A(i,1),A(i,3)) = 1;
%     end 
%     opts={};
%     opts.Nround = 30; %number of EM iterations
%     opts.A_true=A_true;
%     opts.l_true=l_true;
%     opts.tol = 10^-6;
%     disp('Running CNMF-SPA-EM Algorithm....');
%     [A_est_spa_em,l_est_spa_em,timestamps] = run_EM_mod(A_est,l_est,Z,opts);
%     time_em(s,ii)=sum(timestamps);
%     %Resolve permutation ambiguity
%     [rec_error3(s,ii),A_est_p,l_est_p] = getMSE_nomre(A_est_spa_em,l_est_spa_em,A_true,l_true);
%     disp(['Mean Estimation Error - CNMF-SPA-EM: ',num2str(rec_error3(s,ii))]);
%     %disp(['Mean Relative Error - CNMF-SPA-EM: ',num2str(mre_error3(s,ii))]);
%     
%     
%     %% RAND-EM
%     %Run EM algorithm
%     [A] = convert_for_comp(dataset');
%     Z = cell(N,1);
%     for n=1:N
%         Z{n}=zeros(n_items,I(n));
%     end
%     for i = 1:size(A,1)
%         Z{A(i,2)}(A(i,1),A(i,3)) = 1;
%     end
%     opts={};
%     opts.Nround = 30; %number of EM iterations
%     opts.A_true=A_true;
%     opts.l_true=l_true;
%     opts.tol = 10^-6;
%     %l0=1/F_true*ones(F_true,1);
%     disp('Running RAND-EM Algorithm....');
%     [A_est_rand_em,l_est_rand_em,timestamps] = run_EM_mod(A0,l0,Z,opts);
%     time_rand_em(s,ii)=sum(timestamps);
%     %Resolve permutation ambiguity
%     [rec_error4(s,ii),A_est_p,l_est_p] = getMSE_nomre(A_est_rand_em,l_est_rand_em,A_true,l_true);
%     disp(['Mean Estimation Error - RAND-EM: ',num2str(rec_error4(s,ii))]);
%     %disp(['Mean Relative Error - RAND-EM: ',num2str(mre_error4(s,ii))]);
% 
%     %% Tensor Decomposition
%     %Run the LSTensor algorithm
%     [A0,l0]= gen_PMF_factors(I,F_true);
%     %l0=1/F_true*ones(F_true,1);
%     % Algorithm options
%     opts2 = {};
%     for n = 1:N
%         opts2.constraint{n} = 'simplex_col';
%     end
%     opts2.constraint{N+1}   = 'simplex';
%     opts2.max_iter = 200;
%     opts2.rho = [10 10]; 
%     opts2.A0=A0;opts2.l0=l0;
%     opts2.marg    = marg_t;
%     opts2.computeCost=1;
%     opts2.computeMSE=0;
%     opts2.computeCostInterval=20;
%     opts2.tol_impr = 1e-6;
%    disp('Running CTD Algorithm....');
%     [A_est_LS,l_est_LS,Out_LS] = N_CTF_AO_ADMM_Paper_Notation(Y_t,I,F_true,opts2);
%     temp=cumsum(Out_LS.time_instants);
%     time_t(s,ii)=temp(end);
%     %Resolve permutation ambiguity
%    [rec_error5(s,ii),A_est_LS_p,l_est_LS_p] = getMSE_nomre(A_est_LS,l_est_LS,A_true,l_true);
%    disp(['Mean Estimation Error - CTD: ',num2str(rec_error5(s,ii))]);
%    %disp(['Mean Relative Error - CTD: ',num2str(mre_error5(s,ii))]);
%    
%    %% TensorDecomp and EM
%     opts={};
%     opts.Nround = 30; %number of EM iterations
%     opts.A_true=A_true;
%     opts.l_true=l_true;
%     opts.tol = 10^-6;
%     disp('Running CTD-EM Algorithm....');
%     [A_est,l_est,timestamps] = run_EM_mod(A_est_LS,l_est_LS,Z,opts);
%     time_em_ctd(s,ii)=sum(timestamps);
%     %Resolve permutation ambiguity
%     [rec_error6(s,ii),A_est_p,l_est_p] = getMSE_nomre(A_est,l_est,A_true,l_true);
%     disp(['Mean Estimation Error - CTD-EM: ',num2str(rec_error6(s,ii))]);
%     %disp(['Mean Relative Error - CTD-EM: ',num2str(mre_error6(s,ii))]);
% 
    %% Yeredor ML AO
   opts={};
   opts.marg= marg_t;
   opts.max_iter=20;
   opts.tol=10^-6;
   [A_est_mle_ad,l_est_mle_ad,Out_AD] = run_sweeps(Y_t,A0,l0,opts);
  [rec_error7(s,ii),A_est_p,l_est_p] = getMSE_nomre(A_est_mle_ad,l_est_mle_ad,A_true,l_true);
    temp=cumsum(Out_AD.time_instants);
    time_mle_ad(s,ii)=temp(end);
   disp(['Mean Estimation Error - MLE-AD: ',num2str(rec_error7(s,ii))]);
   %disp(['Mean Relative Error - MLE-AD: ',num2str(mre_error7(s,ii))]);
   
   %% Yeredor ML AO
    opts={};
    opts.Nround = 30; %number of EM iterations
    opts.A_true=A_true;
    opts.l_true=l_true;
    opts.tol = 10^-6;
    [A_est_mle_em,l_est_mle_em,timestamps] = run_EM_mod(A_est_mle_ad,l_est_mle_ad,Z,opts);
    time_mle_em(s,ii)=sum(timestamps);
    %Resolve permutation ambiguity
    [rec_error8(s,ii),mre_error8(s,ii),A_est_p,l_est_p] = getMSE_nomre(A_est_mle_em,l_est_mle_em,A_true,l_true);
    disp(['Mean Estimation Error - MLE-EM: ',num2str(rec_error8(s,ii))]);
    disp(['Mean Relative Error - MLE-EM: ',num2str(mre_error8(s,ii))]);
   

   %A_ml = computeML(dataset,F_true,I);   
%    A_ml = computeML(dataset_ML,F_true,I);
%    [rec_error7(s,ii),mre_error7(s,ii),A_est_LS_p,l_est_LS_p] = getMSE(A_ml,l_true,A_true,l_true,1);
%    disp(['Mean Estimation Error - ML: ',num2str(rec_error7(s,ii))]);
%    disp(['Mean Relative Error - ML: ',num2str(mre_error7(s,ii))]);

end
ii=ii+1;
end
end
Nitems = {'N=1000';'N=10000';'N=100000';'N=1000000';'time(s)'};
 %Nitems = {'N=1000';'time(s)'};
 CNMF_SPA = [nanmean(rec_error1,1)';mean(time_mspa,"all")];
 CNMF_OPT = [nanmean(rec_error2,1)';mean(time_m,"all")];
 SPA_EM = [nanmean(rec_error3,1)';mean(time_em,"all")];
 RAND_EM = [nanmean(rec_error4,1)';mean(time_rand_em,"all")];
 CTD = [nanmean(rec_error5,1)';mean(time_t,"all")];
 CTD_EM = [nanmean(rec_error6,1)';mean(time_em_ctd,"all")];
 %MLE_AD = [nanmean(rec_error7,1)';mean(time_mle_ad,"all")];
 %MLE_EM = [nanmean(rec_error8,1)';mean(time_mle_em,"all")];

 T = table(CNMF_SPA,CNMF_OPT,SPA_EM,RAND_EM,CTD,CTD_EM,...
    'RowNames',Nitems);
writetable(T,'table_N_15_I_10_F_10_M_3_p50.csv','WriteRowNames',true)

% Nitems = {'N=1000';'N=10000';'N=100000';'N=1000000';'time(s)'};
%  %Nitems = {'N=1000';'time(s)'};
%  CNMF_SPA = [nanmean(mre_error1,1)';mean(time_mspa,"all")];
%   CNMF_OPT = [nanmean(mre_error2,1)';mean(time_m,"all")];
%  SPA_EM = [nanmean(mre_error3,1)';mean(time_em,"all")];
%  RAND_EM = [nanmean(mre_error4,1)';mean(time_rand_em,"all")];
%  CTD = [nanmean(mre_error5,1)';mean(time_t,"all")];
%  CTD_EM = [nanmean(mre_error6,1)';mean(time_em_ctd,"all")];
%  MLE_AD = [nanmean(mre_error7,1)';mean(time_mle_ad,"all")];
%  MLE_EM = [nanmean(mre_error8,1)';mean(time_mle_em,"all")];
%  T = table(CNMF_SPA,CNMF_OPT,SPA_EM,RAND_EM,CTD,CTD_EM,MLE_AD,MLE_EM,...
%     'RowNames',Nitems);
% writetable(T,'table_mre_N_5_I_10_F_7_M_3_p50_eps_0_3.csv','WriteRowNames',true)

save data_N_15_I_10_F_10_M_3_p50
%     figure(3)
%    loglog(n_items_list,nanmean(rec_error1,1),'-or','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(rec_error2,1),'-*b','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(rec_error3,1),'-dc','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(rec_error4,1),'-<m','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(rec_error5,1),'->g','LineWidth',2)







  
%     figure(3)
%    loglog(n_items_list,nanmean(rec_error1,1),'-or','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(rec_error2,1),'-*b','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(rec_error3,1),'-dc','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(rec_error4,1),'-<m','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(rec_error5,1),'->g','LineWidth',2)


%     grid on;
%         title(['N= ',num2str(N),', I=',num2str(I(1)),', F=',num2str(F_true),', M=',num2str(M)]);
%     xlabel("Number of items");
%     ylabel("MSE");
%     legend("CNMF-SPA","CNMF-OPT","SPA-EM","RAND-EM","CTD");
%     set(gca,'fontsize',14);
%     
%     
%     figure(4)
%    loglog(n_items_list,nanmean(mre_error1,1),'-or','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(mre_error2,1),'-*b','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(mre_error3,1),'-dc','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(mre_error4,1),'-<m','LineWidth',2)
%     hold on;
%     loglog(n_items_list,nanmean(mre_error5,1),'->g','LineWidth',2)
% 
% 
%     grid on;
%         title(['N= ',num2str(N),', I=',num2str(I(1)),', F=',num2str(F_true),', M=',num2str(M)]);
%     xlabel("Number of items");
%     ylabel("MRE");
%     legend("CNMF-SPA","CNMF-OPT","SPA-EM","RAND-EM","CTD");
%     set(gca,'fontsize',14);
    
    
