clc
clearvars
close all

addpath functions
addpath Datasets
addpath functions/tensorlab

load User_Movies
Ratings = round(Ratings); % !!! I have rounded it

MIN_NUMBER_MOVIES       = 5;
PERC_DATA_HIDDEN_test   = 0.3; % hide PERC_DATA_HIDDEN_test of the data for test
PERC_DATA_HIDDEN_valid  = 0.2; % hide PERC_DATA_HIDDEN_valid of the data for validation

 %Movie_ids_most_popular_action_20 = [110 165 260 377 380 457 589 592 648
 %780 1196 1198 1210 1240 1580 2571 3578 4993 5952 7153 2153 1036 3805 2402
 %2823 2818 20 9 71];2115
 Movie_ids_most_popular_action_30 = [110 165 260 377 380 457 589 592 648 780 1196 1198 1210 1240 1580 2571 3578 4993 5952 7153 2153 1036 3805 2402 2823 2818 1370 1036 1291 2115];
 Movie_ids_most_popular_animat_30 = [1 364 551 588 594 595 2355 2700 2987 3114 3751 4306 4886 5218 5618 6377 8360 8961 60069 68954 13 48 558 596 631 720 1022 1029 1032 1033];
 Movie_ids_most_popular_romanc_30 = [39 253 356 357 367 539 587 588 595 597 912 1197 1265 1721 1923 2706 3996 4973 4995 7361 1088 3824 3829 920 2671 339 215 1057 499 605];

Movie_ids = Movie_ids_most_popular_animat_30;

M_rating   = cell(length(Movie_ids),1);
M_r_sparse = cell(length(Movie_ids),1);

for id = 1: length(Movie_ids)
    M_rating{id} = Ratings(Ratings(:,2)== Movie_ids(id),:); M_rating{id}(:,2) = 1;
end
max_elements = 0;
for id = 1: length(Movie_ids)
    M_r_sparse{id} = full(spconvert(M_rating{id}));
    max_elements = max(max_elements,size(M_r_sparse{id},1));
end
for id = 1: length(Movie_ids)
    M_r_sparse{id} = [M_r_sparse{id}; zeros(max_elements - size(M_r_sparse{id},1),1)];
end
User_Movies = zeros(max_elements,length(Movie_ids));
for id = 1: length(Movie_ids)
    User_Movies(:,id) = M_r_sparse{id};
end
Dataset = User_Movies(sum(User_Movies~=0,2)>=MIN_NUMBER_MOVIES,:);
nnz(Dataset)/numel(Dataset);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_sim = 20;
I     = ones(1,size(Dataset,2))*5;
N     = length(I);
F_spa_list    = [5,10,15,20,25];
F_tensor_list=[10,15,20,25];
EM_iteration_list=[5,10,15,20,25,30];


error1 = zeros(n_sim,1);
error2 = zeros(n_sim,1);
error3 = zeros(n_sim,1);
error4 = zeros(n_sim,1);

lbd     =  [10];
BMF_F   =  [5,10,15,20];

rmse = zeros(n_sim,1);
mae  = zeros(n_sim,1);

my_rmse = zeros(n_sim,1);
my_mae  = zeros(n_sim,1);

RMSE_global = zeros(n_sim,1);
RMSE_user = zeros(n_sim,1);
RMSE_movie = zeros(n_sim,1);
MAE_global = zeros(n_sim,1);
MAE_user = zeros(n_sim,1);
MAE_movie = zeros(n_sim,1);

%RMSE_BMF = zeros(n_sim,length(BMF_F));
%MAE_BMF  = zeros(n_sim,length(BMF_F));

% er_rmse = zeros(length(marg),length(F),n_sim);
% er_mae = zeros(length(marg),length(F),n_sim);

for sim=1:n_sim
    nz_ind = find(Dataset>0);
    hidden_elems = datasample(1:length(nz_ind),round((PERC_DATA_HIDDEN_test + PERC_DATA_HIDDEN_valid)*length(nz_ind)),'Replace',false);
    nnz_hidden = nz_ind(hidden_elems);
    
    nnz_hidden_test  = nnz_hidden(1:round( PERC_DATA_HIDDEN_test/(PERC_DATA_HIDDEN_test+PERC_DATA_HIDDEN_valid)*length(nnz_hidden)));
    nnz_hidden_valid = nnz_hidden(length(nnz_hidden_test)+1:end);
    
    TrainSet = Dataset;
    TrainSet(nnz_hidden) = 0;
    
    fprintf('----Dataset Statistics for training set ---- \n');
    fprintf('Number of users    : %d \n', size(TrainSet,1));
    fprintf('Number of movies   : %d \n', size(TrainSet,2));
    fprintf('Number of ratings  : %d \n', nnz(TrainSet));
    fprintf('Density            : %f \n', nnz(TrainSet)/numel(TrainSet));
    fprintf('Size of training set    : %f \n', nnz(TrainSet));
    fprintf('Size of validation set  : %f \n', length(nnz_hidden_valid));
    fprintf('Size of test set        : %f \n', length(nnz_hidden_test));
    
    [r1,c1] = find(TrainSet == 0);       % training set indices
    for i = 1: length(nnz_hidden_valid)
        [r3(i),c3(i)] = ind2sub(size(Dataset),nnz_hidden_valid(i));  % validation set indices
        P3(i) = Dataset(r3(i),c3(i));                                % validation set values
    end
    for i = 1: length(nnz_hidden_test)
        [r4(i),c4(i)] = ind2sub(size(Dataset),nnz_hidden_test(i)); % test set indices
        P4(i) = Dataset(r4(i),c4(i));                              % test set values
    end
% %     %     %%%%%%%%%%%%%%%%%%%%%%% Biased Matrix Factorization %%%%%%%%%%%%%%%
        for l = 1:length(lbd)
            fprintf('Lambda used: %d \n',lbd(l))
            for f = 1: length(BMF_F)
                fprintf('Rank used: %d \n',BMF_F(f))
                a=tic;
                [A,B,cost,RMSE_v,RMSE_t,MAE_t] = BiasedMatrixFactorization(TrainSet,r1,c1,P3,r3,c3,P4,r4,c4,lbd(l),BMF_F(f));
                b=toc(a);
                time_bmf(sim,f)=b;
                RMSE_BMF_v(f) = RMSE_v;
                RMSE_BMF_t(f) = RMSE_t;
                MAE_BMF_t(f) = MAE_t;
            end
        end
        [~,ind]=min(RMSE_BMF_v);
        RMSE_BMF(sim)=RMSE_BMF_t(ind);
        MAE_BMF(sim)=MAE_BMF_t(ind);
%         
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%% Simple Methods %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Global Average, User Average, Item Average
    global_average = sum(TrainSet(TrainSet~=0))/nnz(TrainSet);
    user_average   = sum(TrainSet,2)./sum(TrainSet~=0,2);
    movie_average  = sum(TrainSet,1)./sum(TrainSet~=0,1);
    
    count      = 0;
    Err_global = 0;
    Err_user   = 0;
    Err_movie  = 0;
    
    ErrMAE_global = 0;
    ErrMAE_user   = 0;
    ErrMAE_movie  = 0;
    
    for i = 1:length(P4)
        if any(TrainSet(r4(i),:)~=0)
            Err_global = Err_global + (global_average - P4(i))^2;
            Err_user   = Err_user   + (user_average(r4(i)) - P4(i))^2;
            Err_movie  = Err_movie  + (movie_average(c4(i)) - P4(i))^2;
            
            ErrMAE_global = ErrMAE_global + abs(global_average - P4(i));
            ErrMAE_user   = ErrMAE_user   + abs(user_average(r4(i)) - P4(i));
            ErrMAE_movie  = ErrMAE_movie  + abs(movie_average(c4(i)) - P4(i));
            count = count + 1;
        end
    end
    RMSE_global(sim) = sqrt(Err_global/count)
    RMSE_user(sim)   = sqrt(Err_user/count)
    RMSE_movie(sim)  = sqrt(Err_movie/count)
    MAE_global(sim)  = ErrMAE_global/count
    MAE_user(sim)    = ErrMAE_user/count
    MAE_movie(sim)   = ErrMAE_movie/count
    fprintf('Number of ratings predicted  : %d instead of %d \n', count, length(nnz_hidden_test));
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% SPA-Our Approach %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fprintf('Marginal used: 2')
    marg  = combnk(1:length(I),2);
    marg = num2cell(marg,2);
    Y = get_obs_marg_fast(TrainSet,marg,I);
    A_LS_spa = {};
    l_LS_spa ={};
    rmse_valid_spa=[];
    mae_valid_spa=[];
    for f= 1:length(F_spa_list)
        R = F_spa_list(f);
        fprintf('Rank used: %d \n',R)
        %if(R > min(I))
        max_iter=300;
        tol=10^-6;
        M=5;
        opts={};
        opts.marg=marg;
        opts.max_iter=max_iter;
        opts.tol=tol;
        disp('Running SingleSPA Algorithm....');
        a=tic;
        [A_est,l_est] = SingleSPA(Y,I,M,R,opts);
        b=toc(a);
        time_spa(sim,f)=b;
        for i=1:N
           A_est{i}(isnan(A_est{i}))=0;
        end
            
        [rmse_valid_spa(f),mae_valid_spa(f),pr_rating,num_nan] = complete_rating_matrix(TrainSet,P3,r3,c3,A_est,l_est);


        A_LS_spa{f}=A_est;
        l_LS_spa{f}=l_est;
    end
    [~,ind]=min(rmse_valid_spa);
    F_sel=F_spa_list(ind);
    A_sel=A_LS_spa{ind};
    l_sel=l_LS_spa{ind};
    [rmse_spa(sim),mae_spa(sim),pr_rating,num_nan] = complete_rating_matrix(TrainSet,P4,r4,c4,A_sel,l_sel);

%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Matrix KL Method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      A_KL_m = {};
%      l_KL_m ={};
%      rmse_valid_km=[];
%      mae_valid_km=[];
%      for f= 1:length(F_spa_list)
%         R = F_spa_list(f);
%         fprintf('Rank used: %d \n',R)
%         % Algorithm options
%         opts2 = {};
%         for n = 1:N
%             opts2.constraint{n} = 'simplex_col';
%         end
%         opts2.constraint{N+1}   = 'simplex';
%         opts2.max_iter = 200;
%         opts2.rho = [10 10]; 
%         A_init=A_LS_spa{f};
%         for i=1:length(A_init)
%            G =A_init{i};
%            G = max( G, 10^-6 );
%            t=sum(G,1);
%            G = G*diag(1./t);
%            A_init{i}=G;
%         end
%         opts2.A0=A_init;opts2.l0=1/R*ones(R,1);%l_LS_spa{f};
%         opts2.marg    = marg;
%         opts2.computeCost=1;
%         opts2.computeMSE=0;
%         opts2.computeCostInterval=20;
%         opts2.tol_impr = 1e-4;
%         disp('Running KL_Matrix Algorithm....');
%         [A_est_km,l_est_km,Out_km] = N_CTF_AO_KL(Y,I,R,opts2);
%         temp=cumsum(Out_km.time_instants);
%         time_km(sim,f)=temp(end);
%         [rmse_valid_km(f),mae_valid_km(f),pr_rating,num_nan] = complete_rating_matrix(TrainSet,P3,r3,c3,A_est_km,l_est_km);
%         A_KL_m{f}=A_est_km;
%         l_KL_m{f}=l_est_km;
%     end
%     [~,ind]=min(rmse_valid_km);
%     F_sel=F_spa_list(ind);
%     A_sel=A_KL_m{ind};
%     l_sel=l_KL_m{ind};
%     [rmse_km(sim),mae_km(sim),pr_rating,num_nan] = complete_rating_matrix(TrainSet,P4,r4,c4,A_sel,l_sel);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%SPA-EM Method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     A_em = {};
     l_em ={};
     rmse_valid_em=[];
     mae_valid_em=[];
     for f= 1:length(F_spa_list)
        R = F_spa_list(f);
        fprintf('Rank used: %d \n',R)
        iter=1;
        for Nround=EM_iteration_list
            % Algorithm 
            [A] = convert_for_comp(TrainSet');
            trainset_len= size(TrainSet,1);
            Z = cell(N,1);
            for n=1:N
                Z{n}=zeros(trainset_len,I(n));
            end
            for i = 1:size(A,1)
                Z{A(i,2)}(A(i,1),A(i,3)) = 1;
            end
            disp('Running EM Algorithm....');
            a=tic;
            %Nround = 10; %number of EM iterations
            [A_est_em,l_est_em] = run_EM(A_LS_spa{f},l_LS_spa{f},Z,Nround);
            b=toc(a);
            time_em(sim,f)=b;
            [rmse_valid_em(iter,f),mae_valid_em(iter,f),pr_rating,num_nan] = complete_rating_matrix(TrainSet,P3,r3,c3,A_est_em,l_est_em);
            A_em{iter,f}=A_est_em;
            l_em{iter,f}=l_est_em;
            iter=iter+1;
        end        
     end
    max_val=max(rmse_valid_em,[],'all');
    [ind1,ind2]=find(rmse_valid_em==max_val);
    F_sel=F_spa_list(ind2(1));
%     A_sel=A_em{ind1(1),ind2(1)};
%     l_sel=l_em{ind1(1),ind2(1)};
     [A_sel,l_sel] = run_EM(A_LS_spa{ind2(1)},l_LS_spa{ind2(1)},Z,Nround);
    [rmse_em(sim),mae_em(sim),pr_rating,num_nan] = complete_rating_matrix(TrainSet,P4,r4,c4,A_sel,l_sel);


    %%%%%%%%%%%%%%%%%%%%Tesnor Method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Marginal used: 3')
    opts={};
    opts.marg  = combnk(1:length(I),3);
    opts.marg = num2cell(opts.marg,2);
    Y = get_obs_marg_fast(TrainSet,opts.marg,I);
    A_LS_t = {};
    l_LS_t ={};
    rmse_valid_t=[];
    mae_valid_t=[];
    for f= 1:length(F_tensor_list)
        R = F_tensor_list(f);
        fprintf('Rank used: %d \n',R)
        [opts.A0,l0]   = gen_PMF_factors(I,R);
        opts.l0  = 1/R*ones(R,1);
        %opts.A0=A_est;
        %opts.l0=l_est;
        opts.max_iter = 200;
        opts.tol_impr = 1e-4;
        for n = 1:N
            opts.constraint{n} = 'simplex_col';
        end
        opts.constraint{N+1}   = 'simplex';

        opts.TrainSet = TrainSet;
        opts.P3 = P3;
        opts.r3 = r3;
        opts.c3 = c3;

        opts.P4 = P4;
        opts.r4 = r4;
        opts.c4 = c4;

        opts.rho = [10 10];
        [A,lambda,Out] = N_CTF_AO_ADMM_Paper_Notation(Y,I,R,opts);
        temp=cumsum(Out.time_instants);
        time_t(sim,f)=temp(end);
        [rmse_valid_t(f),mae_valid_t(f),pr_rating,num_nan] = complete_rating_matrix(TrainSet,P3,r3,c3,A,lambda);
        A_LS_t{f}=A;
        l_LS_t{f}=lambda;
    end
    [~,ind]=max(rmse_valid_t);
    F_sel=F_tensor_list(ind);
    A_sel=A_LS_t{ind};
    l_sel=l_LS_t{ind};
    [rmse_t(sim),mae_t(sim),pr_rating,num_nan] = complete_rating_matrix(TrainSet,P4,r4,c4,A_sel,l_sel);
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CTD-EM Method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     A_em_ctd = {};
     l_em_ctd ={};
     rmse_valid_em_ctd=[];
     mae_valid_em_ctd=[];
     for f= 1:length(F_tensor_list)
        R = F_tensor_list(f);
        fprintf('Rank used: %d \n',R)
        iter=1;
        for Nround=EM_iteration_list
            % Algorithm 
            [A] = convert_for_comp(TrainSet');
            trainset_len= size(TrainSet,1);
            Z = cell(N,1);
            for n=1:N
                Z{n}=zeros(trainset_len,I(n));
            end
            for i = 1:size(A,1)
                Z{A(i,2)}(A(i,1),A(i,3)) = 1;
            end
            disp('Running EM Algorithm....');
            a=tic;
            %Nround = 10; %number of EM iterations
            [A_est_em_ctd,l_est_em_ctd] = run_EM(A_LS_t{f},l_LS_t{f},Z,Nround);
            b=toc(a);
            time_em_ctd(sim,f)=b;
            [rmse_valid_em_ctd(iter,f),mae_valid_em_ctd(iter,f),pr_rating,num_nan] = complete_rating_matrix(TrainSet,P3,r3,c3,A_est_em_ctd,l_est_em_ctd);
            A_em_ctd{iter,f}=A_est_em_ctd;
            l_em_ctd{iter,f}=l_est_em_ctd;
            iter=iter+1;
        end        
     end
    max_val=max(rmse_valid_em_ctd,[],'all');
    [ind1,ind2]=find(rmse_valid_em_ctd==max_val);
    F_sel=F_tensor_list(ind2(1));
%     A_sel=A_em{ind1(1),ind2(1)};
%     l_sel=l_em{ind1(1),ind2(1)};
     [A_sel,l_sel] = run_EM(A_LS_t{ind2(1)},l_LS_t{ind2(1)},Z,Nround);
    [rmse_em_ctd(sim),mae_em_ctd(sim),pr_rating,num_nan] = complete_rating_matrix(TrainSet,P4,r4,c4,A_sel,l_sel);

            
end
rows = {'RMSE';'MAE';'STD-RMSE';'STD-MAE';'time(s)'};
 CNMF_SPA = [nanmean(rmse_spa);nanmean(mae_spa);std(rmse_spa);std(mae_spa);mean(time_spa,"all")];
 %CNMF_OPT = [nanmean(rmse_km);nanmean(mae_km);std(rmse_km);std(mae_km);mean(time_km,"all")];
 SPA_EM = [nanmean(rmse_em);nanmean(mae_em);std(rmse_em);std(mae_em);mean(time_em,"all")];
 CTD = [nanmean(rmse_t);nanmean(mae_t);std(rmse_t);std(mae_t);mean(time_t,"all")];
 CTD_EM = [nanmean(rmse_em_ctd);nanmean(mae_em_ctd);std(rmse_em_ctd);std(mae_em_ctd);mean(time_em_ctd,"all")];
 BMF = [nanmean(RMSE_BMF);nanmean(MAE_BMF);std(RMSE_BMF);std(MAE_BMF);mean(time_bmf,"all")];
 Global_average = [nanmean(RMSE_global);nanmean(MAE_global);std(RMSE_global);std(MAE_global);mean(time_bmf,"all")];
 User_average   = [nanmean(RMSE_user);nanmean(MAE_user);std(RMSE_user);std(MAE_user);mean(time_bmf,"all")];
 Movie_average  = [nanmean(RMSE_movie);nanmean(MAE_movie);std(RMSE_movie);std(MAE_movie);mean(time_bmf,"all")];

 T = table(CNMF_SPA,SPA_EM,CTD,CTD_EM,BMF,Global_average,User_average,Movie_average,...
    'RowNames',rows);
writetable(T,'table_animation.csv','WriteRowNames',true)
save data_animation
% disp('########Results##############');
% disp(['RMSE - SPA ',num2str(mean(rmse_spa))]);
% disp(['RMSE - KLMatrix ',num2str(mean(rmse_km))]);
% disp(['RMSE - SPA-EM ',num2str(mean(rmse_em))]);
% disp(['RMSE - LSTensor ',num2str(mean(rmse_t))]);
% disp(['RMSE - BMF ',num2str(mean(RMSE_BMF))]);
% disp(['RMSE - Global Avergae ',num2str(mean(RMSE_global))]);
% disp(['RMSE - User Avergae ',num2str(mean(RMSE_user))]);
% disp(['RMSE - Movie Avergae ',num2str(mean(RMSE_movie))]);
% disp('########Results##############');
% disp(['MAE - SPA ',num2str(mean(mae_spa))]);
% disp(['MAE - KLMatrix ',num2str(mean(mae_km))]);
% disp(['MAE - SPA-EM ',num2str(mean(mae_em))]);
% disp(['MAE - LSTensor ',num2str(mean(mae_t))]);
% disp(['MAE - BMF ',num2str(mean(MAE_BMF))]);
% disp(['MAE - Global Avergae ',num2str(mean(MAE_global))]);
% disp(['MAE - User Avergae ',num2str(mean(MAE_user))]);
% disp(['MAE - Movie Avergae ',num2str(mean(MAE_movie))]);
% disp('########Results##############');
% disp(['TIME - SPA ',num2str(mean(time_spa,"all"))]);
% disp(['TIME - KLMatrix ',num2str(mean(time_km,"all"))]);
% disp(['TIME - SPA-EM ',num2str(mean(time_em,"all"))]);
% disp(['TIME - LSTensor ',num2str(mean(time_t,"all"))]);
% disp(['TIME - BMF ',num2str(mean(time_bmf,"all"))]);
% 
% disp('########Results##############');
% disp(['RMSE-STD - SPA ',num2str(std(rmse_spa))]);
% disp(['RMSE-STD - KLMatrix ',num2str(std(rmse_km))]);
% disp(['RMSE-STD - SPA-EM ',num2str(std(rmse_em))]);
% disp(['RMSE-STD - LSTensor ',num2str(std(rmse_t))]);
% disp(['RMSE-STD - BMF ',num2str(std(RMSE_BMF))]);
% disp(['RMSE-STD - Global Avergae ',num2str(std(RMSE_global))]);
% disp(['RMSE-STD - User Avergae ',num2str(std(RMSE_user))]);
% disp(['RMSE-STD - Movie Avergae ',num2str(std(RMSE_movie))]);
% 
% disp('########Results##############');
% disp(['MAE-STD - SPA ',num2str(std(mae_spa))]);
% disp(['MAE-STD - KLMatrix ',num2str(std(mae_km))]);
% disp(['MAE-STD - SPA-EM ',num2str(std(mae_em))]);
% disp(['MAE-STD - LSTensor ',num2str(std(mae_t))]);
% disp(['MAE-STD - BMF ',num2str(std(MAE_BMF))]);
% disp(['MAE-STD - Global Avergae ',num2str(std(MAE_global))]);
% disp(['MAE-STD - User Avergae ',num2str(std(MAE_user))]);
% disp(['MAE-STD - Movie Avergae ',num2str(std(MAE_movie))]);
% writetable(T,'table_action.csv','WriteRowNames',true)
% save data_action