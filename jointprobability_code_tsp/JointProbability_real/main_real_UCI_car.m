clc
clearvars
close all

addpath functions
addpath functions/tensorlab
addpath Datasets
addpath functions_mle

%Test Variables
[Data,I]    = LoadDataset_car();
[rows, col] = size(Data);
N           = length(I);
F_list_SPA=[2,3,4,5];
F_list_SPA1=[2,3,4,5];
F_list_Tensor=[2,3,4,5];
EM_iteration_list=[5,10,15,20,25,30];
M=5;

% Number of simulations
n_sim = 10;


%GD parameters
max_iter=300;
tol=10^-6;

%Algorithms
SINGLE_SPA=1;
SPA_EM=1;
LS_TENSOR=1;
KL_MATRIX=1;
SVM=1;
LOG_REGRESSION=1;
LINEAR_REG=1;
NEURAL_NET=1;
NAIVE_BAYES=1;
SVM_RBF=1;
CTD_EM=0;
MLE_AD=1;
MLE_EM=1;

no_of_class = max(Data(:,1));

for s = 1:n_sim

    [trainIndices,validIndices,testIndices] = dividerand(rows,0.5,0.3,0.2);

    Data_train   = Data(trainIndices,2:end);
    Data_test    = Data(testIndices,2:end);
    Data_valid   = Data(validIndices,2:end);


    Labels_train = Data(trainIndices,1);
    Labels_test  = Data(testIndices,1);
    Labels_valid = Data(validIndices,1);

    TrainSet      = [Labels_train Data_train];
    TestSet       = [Labels_test  Data_test];
    ValidSet      = [Labels_valid Data_valid];

    a=tic;
    marg    = num2cell(combnk(1:length(I),2),2);
    Y = get_obs_marg_fast(TrainSet,marg,I);        % Get the lower-dimensional tensors
    b=toc(a);
    time_marg2(s) = b;
    
    a=tic;
    marg_t    = num2cell(combnk(1:length(I),3),2);
    Y_t = get_obs_marg_fast(TrainSet,marg_t,I);        % Get the lower-dimensional tensors
    b=toc(a);
    time_marg3(s)=b;
    
    len=length(Labels_test);
    len_v=length(Labels_valid);
    if(SINGLE_SPA)
        fi=1;
        A_SPA = {};
        l_SPA ={};
        accuracy_valid_spa=[];
        for F_true = F_list_SPA
            %Run the initialization using singleSPA algorithm 
            opts={};
            opts.marg=marg;
            opts.max_iter=max_iter;
            opts.tol=tol;
            disp('Running SingleSPA Algorithm....');
            a=tic;
            [A_est,l_est] = SingleSPA(Y,I,M,F_true,opts);
            b=toc(a);

            time_spa(s,fi)=b;
            for i=1:N
               A_est{i}(isnan(A_est{i}))=0;
            end
            for i=1:length(A_est)
               G =A_est{i};
               G = max( G, 10^-6 );
               t=sum(G,1);
               G = G*diag(1./t);
               A_est{i}=G;
            end
            [Labels,ct]=predict_parafac_multiclass(A_est,l_est,Data_valid,I,F_true);
            accuracy_valid_spa(fi)=sum(Labels==Labels_valid)/len_v
            A_SPA{fi}=A_est;
            l_SPA{fi}=l_est;
            fi=fi+1;
        end
        [~,ind]=max(accuracy_valid_spa);
        F_sel=F_list_SPA(ind);
        A_sel=A_SPA{ind};
        l_sel=l_SPA{ind};
        [Labels,ct]=predict_parafac_multiclass(A_sel,l_sel,Data_test,I,F_sel);
        accuracy_test_spa(s)=sum(Labels==Labels_test)/len
    end
    
    if(SPA_EM)
        fi=1;
        A_em = {};
        l_em ={};
        accuracy_valid_em=[];
        for F_true = F_list_SPA
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
                %Nround = 20; %number of EM iterations
                [A_est_em,l_est_em] = run_EM(A_SPA{fi},l_SPA{fi},Z,Nround);
                b=toc(a);
                time_em(s,fi)=b;
    %             for i=1:N
    %                A_est_m{i}(isnan(A_est_m{i}))=0;
    %             end
                [Labels_m,ct]=predict_parafac_multiclass(A_est_em,l_est_em,Data_valid,I,F_true);
                accuracy_valid_em(iter,fi)=sum(Labels_m==Labels_valid)/len_v
                A_em{iter,fi}=A_est_em;
                l_em{iter,fi}=l_est_em;
                iter=iter+1;
            end
            fi=fi+1;
        end
        %[~,ind]=max(accuracy_valid_em);
        max_val=max(accuracy_valid_em,[],'all');
        [ind1,ind2]=find(accuracy_valid_em==max_val);
        F_sel=F_list_SPA(ind2(1));
%         A_sel=A_em{ind1(1),ind2(1)};
%         l_sel=l_em{ind1(1),ind2(1)};
        [A_sel,l_sel] = run_EM(A_SPA{ind2(1)},l_SPA{ind2(1)},Z,Nround);

        [Labels,ct]=predict_parafac_multiclass(A_sel,l_sel,Data_test,I,F_sel);
        accuracy_test_em(s)=sum(Labels==Labels_test)/len
    end
    
 
    if(KL_MATRIX)
        fi=1;
        A_KL_m = {};
        l_KL_m ={};
        accuracy_valid_km=[];
        for F_true = F_list_SPA
            %Run the LS algorithm
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
            A_est=A_SPA{fi};
            for i=1:length(A_est)
               G =A_est{i};
               G = max( G, 10^-6 );
               t=sum(G,1);
               G = G*diag(1./t);
               A_est{i}=G;
            end
            opts2.A0=A_est;opts2.l0=l_SPA{fi};
            opts2.marg    = marg;
            opts2.computeCost=1;
            opts2.computeMSE=0;
            opts2.computeCostInterval=5;
            opts2.tol_impr = 1e-6;
            disp('Running KL_Matrix Algorithm....');
            [A_est_km,l_est_km,Out_km] = N_CTF_AO_KL(Y,I,F_true,opts2);
            temp=cumsum(Out_km.time_instants);
            time_km(s,fi)=temp(end);
%             for i=1:N
%                A_est_m{i}(isnan(A_est_m{i}))=0;
%             end
            [Labels_m,ct]=predict_parafac_multiclass(A_est_km,l_est_km,Data_valid,I,F_true);
            accuracy_valid_km(fi)=sum(Labels_m==Labels_valid)/len_v
            A_KL_m{fi}=A_est_km;
            l_KL_m{fi}=l_est_km;
            fi=fi+1;
        end
        [~,ind]=max(accuracy_valid_km);
        F_sel=F_list_SPA(ind);
        A_sel=A_KL_m{ind};
        l_sel=l_KL_m{ind};
        [Labels,ct]=predict_parafac_multiclass(A_sel,l_sel,Data_test,I,F_sel);
        accuracy_test_km(s)=sum(Labels==Labels_test)/len
    end
    


    if(LS_TENSOR)
        fi=1;
        A_LS_t = {};
        l_LS_t ={};
        accuracy_valid_t=[];
        for F_true = F_list_Tensor
            %Run the LSTensor algorithm
            [A0,l0]= gen_PMF_factors(I,F_true);
            %l0=1/F_true*ones(F_true,1);
            % Algorithm options
            opts2 = {};
            for n = 1:N
                opts2.constraint{n} = 'simplex_col';
            end
            opts2.constraint{N+1}   = 'simplex';
            opts2.max_iter = 100;
            opts2.rho = [10 10]; 
            opts2.A0=A0;opts2.l0=l0;
            opts2.marg    = marg_t;
            opts2.computeCost=1;
            opts2.computeMSE=0;
            opts2.computeCostInterval=20;
            opts2.tol_impr = 1e-6;
            disp('Running LS_Tensor Algorithm....');
            [A_est_t,l_est_t,Out_t] = N_CTF_AO_ADMM_Paper_Notation(Y_t,I,F_true,opts2);
            temp=cumsum(Out_t.time_instants);
            time_t(s,fi)=temp(end);
            [Labels_t,ct]=predict_parafac_multiclass(A_est_t,l_est_t,Data_valid,I,F_true);
            accuracy_valid_t(fi)=sum(Labels_t==Labels_valid)/len_v
            A_LS_t{fi}=A_est_t;
            l_LS_t{fi}=l_est_t;
            fi=fi+1;
        end
        [~,ind]=max(accuracy_valid_t);
        F_sel=F_list_Tensor(ind);
        A_sel=A_LS_t{ind};
        l_sel=l_LS_t{ind};
        [Labels,ct]=predict_parafac_multiclass(A_sel,l_sel,Data_test,I,F_sel);
        accuracy_test_t(s)=sum(Labels==Labels_test)/len
    end
    if(CTD_EM)
        fi=1;
        A_em = {};
        l_em ={};
        accuracy_valid_em=[];
        for F_true = F_list_Tensor
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
                %Nround = 20; %number of EM iterations
                [A_est_em,l_est_em] = run_EM(A_LS_t{fi},l_LS_t{fi},Z,Nround);
                b=toc(a);
                time_em_ctd(s,fi)=b;
    %             for i=1:N
    %                A_est_m{i}(isnan(A_est_m{i}))=0;
    %             end
                [Labels_m,ct]=predict_parafac_multiclass(A_est_em,l_est_em,Data_valid,I,F_true);
                accuracy_valid_em(iter,fi)=sum(Labels_m==Labels_valid)/len_v
                A_em{iter,fi}=A_est_em;
                l_em{iter,fi}=l_est_em;
                iter=iter+1;
            end
            fi=fi+1;
        end
        %[~,ind]=max(accuracy_valid_em);
        max_val=max(accuracy_valid_em,[],'all');
        [ind1,ind2]=find(accuracy_valid_em==max_val);
        F_sel=F_list_Tensor(ind2(1));
%         A_sel=A_em{ind1(1),ind2(1)};
%         l_sel=l_em{ind1(1),ind2(1)};
        [A_sel,l_sel] = run_EM(A_LS_t{ind2(1)},l_LS_t{ind2(1)},Z,Nround);

        [Labels,ct]=predict_parafac_multiclass(A_sel,l_sel,Data_test,I,F_sel);
        accuracy_test_em_ctd(s)=sum(Labels==Labels_test)/len
    end
    
    if(MLE_AD)
        fi=1;
        A_MLE_t = {};
        l_MLE_t ={};
        accuracy_valid_mle=[];
        for F_true = F_list_Tensor
            %Run the LSTensor algorithm
            [A0,l0]= gen_PMF_factors(I,F_true);
            opts2 = {};
            opts2.marg    = marg_t;
            opts2.max_iter=20;
            opts2.tol=10^-4;
            disp('Running MLE AD Algorithm....');
            [A_est_t,l_est_t,Out] = run_sweeps(Y_t,A0,l0,opts2);
            temp=cumsum(Out.time_instants);
            time_mle(s,fi)=temp(end);
            [Labels_t,ct]=predict_parafac_multiclass(A_est_t,l_est_t,Data_valid,I,F_true);
            accuracy_valid_mle(fi)=sum(Labels_t==Labels_valid)/len_v
            A_MLE_t{fi}=A_est_t;
            l_MLE_t{fi}=l_est_t;
            fi=fi+1;
        end
        [~,ind]=max(accuracy_valid_mle);
        F_sel=F_list_Tensor(ind);
        A_sel=A_MLE_t{ind};
        l_sel=l_MLE_t{ind};
        [Labels,ct]=predict_parafac_multiclass(A_sel,l_sel,Data_test,I,F_sel);
        accuracy_test_mle(s)=sum(Labels==Labels_test)/len
    end
    if(MLE_EM)
        fi=1;
        A_mle_em = {};
        l_mle_em ={};
        accuracy_valid_mle_em=[];
        for F_true = F_list_Tensor
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
                %Nround = 20; %number of EM iterations
                [A_est_em,l_est_em] = run_EM(A_MLE_t{fi},l_MLE_t{fi},Z,Nround);
                b=toc(a);
                time_mle_em(s,fi)=b;
    %             for i=1:N
    %                A_est_m{i}(isnan(A_est_m{i}))=0;
    %             end
                [Labels_m,ct]=predict_parafac_multiclass(A_est_em,l_est_em,Data_valid,I,F_true);
                accuracy_valid_mle_em(iter,fi)=sum(Labels_m==Labels_valid)/len_v
                A_mle_em{iter,fi}=A_est_em;
                l_mle_em{iter,fi}=l_est_em;
                iter=iter+1;
            end
            fi=fi+1;
        end
        %[~,ind]=max(accuracy_valid_em);
        max_val=max(accuracy_valid_mle_em,[],'all');
        [ind1,ind2]=find(accuracy_valid_mle_em==max_val);
        F_sel=F_list_Tensor(ind2(1));
%         A_sel=A_em{ind1(1),ind2(1)};
%         l_sel=l_em{ind1(1),ind2(1)};
        [A_sel,l_sel] = run_EM(A_MLE_t{ind2(1)},l_MLE_t{ind2(1)},Z,Nround);

        [Labels,ct]=predict_parafac_multiclass(A_sel,l_sel,Data_test,I,F_sel);
        accuracy_test_mle_em(s)=sum(Labels==Labels_test)/len
    end


    if(SVM)
        %SVM
        a=tic;
        if(no_of_class==2)
            SVMModel_1       = fitcsvm(Data_train,Labels_train);
        else
            SVMModel_1       = fitcecoc(Data_train,Labels_train);
        end

        SVMModel_1 = crossval(SVMModel_1,'Holdout',0.10);
        b=toc(a);
        time_svm(s)=b;
        Labels_SVM_1     = predict(SVMModel_1.Trained{1},Data_test);
        accuracy_test_svm(s)=sum(Labels_SVM_1==Labels_test)/len
    end
    
%     if(LINEAR_REG)
%         a=tic;
%         LinearReg_Model = fitcdiscr(Data_train,Labels_train);
%         LinearReg_Model = crossval(LinearReg_Model,'Holdout',0.10);
%         b=toc(a);
%         time_linear_reg(s)=b;
%         Labels_LIN_REG = predict(LinearReg_Model.Trained{1},Data_test);
%         accuracy_test_linearreg(s)=sum(round(Labels_LIN_REG)==Labels_test)/len
%     end
    if(LOG_REGRESSION)
        a=tic;
        LinearReg_Model = mnrfit(Data_train,Labels_train);
        b=toc(a);
        time_linear_reg(s)=b;
         LinearReg_Model = mnrval(LinearReg_Model,Data_test);
%         Labels_LIN_REG = predict(LinearReg_Model,Data_test);
        [~,Labels_LIN_REG] = max(LinearReg_Model,[],2);
        accuracy_test_linearreg(s)=sum(round(Labels_LIN_REG)==Labels_test)/len
    end
    
    if(NEURAL_NET)
        a=tic;
        net=feedforwardnet(10);
        net = train(net,Data_train',Labels_train');
        b=toc(a);
        time_net(s)=b;
        Labels_nnet = net(Data_test');
        accuracy_test_nnet(s)=sum(round(Labels_nnet')==Labels_test)/len
    end
    
    if(SVM_RBF)
        a=tic;
        if(no_of_class==2)
            SVMModel_2       = fitcsvm(Data_train,Labels_train,'KernelFunction','rbf');
        else
            SVMModel_2       = fitcecoc(Data_train,Labels_train,'Learners','kernel');            
        end
        %SVMModel_2 = crossval(SVMModel_2,'Holdout',0.10);
        Labels_SVM_2     = predict(SVMModel_2,Data_test);
        b=toc(a);
        time_svm_rbf(s)=b;
        accuracy_test_svmrbf(s)=sum(Labels_SVM_2==Labels_test)/len
    end
    
    if(NAIVE_BAYES)
       a=tic;
       NB =fitcnb(Data_train,Labels_train,'CategoricalPredictors',1:size(Data_train,2)); %fitclinear
       NB= crossval(NB,'Holdout',0.10);
       b=toc(a);
       time_nb(s)=b;
       Labels_NB        = predict(NB.Trained{1},Data_test);
       accuracy_test_nb(s)=sum(Labels_NB==Labels_test)/len
    end
end

disp(['Accuracy-SPA ',num2str(mean(accuracy_test_spa))]);
disp(['Accuracy-SPA-EM ',num2str(mean(accuracy_test_em))]);
disp(['Accuracy-KL_Matrix ',num2str(mean(accuracy_test_km))]);
disp(['Accuracy-LS_Tensor ',num2str(mean(accuracy_test_t))]);
disp(['Accuracy-SVM ',num2str(mean(accuracy_test_svm))]);
disp(['Accuracy-LinearReg ',num2str(mean(accuracy_test_linearreg))]);
disp(['Accuracy-nnet ',num2str(mean(accuracy_test_nnet))]);
disp(['Accuracy-SVM RBF ',num2str(mean(accuracy_test_svmrbf))]);
disp(['Accuracy-Naive Bayes ',num2str(mean(accuracy_test_nb))]);

disp(['std_Accuracy-SPA ',num2str(std(accuracy_test_spa))]);
disp(['std_Accuracy-SPA-EM ',num2str(std(accuracy_test_em))]);
disp(['std_Accuracy-KL_Matrix ',num2str(std(accuracy_test_km))]);
disp(['std_Accuracy-LS_Tensor ',num2str(std(accuracy_test_t))]);
disp(['std_Accuracy-SVM ',num2str(std(accuracy_test_svm))]);
disp(['std_Accuracy-LinearReg ',num2str(std(accuracy_test_linearreg))]);
disp(['std_Accuracy-nnet ',num2str(std(accuracy_test_nnet))]);
disp(['std_Accuracy-SVM RBF ',num2str(std(accuracy_test_svmrbf))]);
disp(['std_Accuracy-Naive Bayes ',num2str(std(accuracy_test_nb))]);

disp(['Time-SPA ',num2str(mean(time_spa,"all"))]);
disp(['Time-SPA-EM ',num2str(mean(time_em,"all"))]);
disp(['Time-KL_Matrix ',num2str(mean(time_km,"all"))]);
disp(['Time-LS_Tensor ',num2str(mean(time_t,"all"))]);
disp(['Time-SVM ',num2str(mean(time_svm))]);
disp(['Time-LinearReg ',num2str(mean(time_linear_reg))]);
disp(['Time-nnet ',num2str(mean(time_net))]);
disp(['Time-SVM RBF ',num2str(mean(time_svm_rbf))]);
disp(['Time-Naive Bayes ',num2str(mean(time_nb))]);

disp(['Time-Marg2 ',num2str(mean(time_marg2))]);
disp(['Time-Marg3 ',num2str(mean(time_marg3))]);

rows = {'Accuracy';'STD_Accuracy';'time(s)'};
 CNMF_SPA = [nanmean(accuracy_test_spa);std(accuracy_test_spa);mean(time_spa,"all")];
 CNMF_OPT = [nanmean(accuracy_test_km);std(accuracy_test_km);mean(time_km,"all")];
 SPA_EM = [nanmean(accuracy_test_em);std(accuracy_test_em);mean(time_em,"all")];
 CTD = [nanmean(accuracy_test_t);std(accuracy_test_t);mean(time_t,"all")];
 %CTD_EM = [nanmean(accuracy_test_em_ctd);std(accuracy_test_em_ctd);mean(time_em_ctd,"all")];
 MLE_AD = [nanmean(accuracy_test_mle);std(accuracy_test_mle);mean(time_mle,"all")];
 MLE_EM = [nanmean(accuracy_test_mle_em);std(accuracy_test_mle_em);mean(time_mle_em,"all")];
 SVM = [nanmean(accuracy_test_svm);std(accuracy_test_svm);mean(time_svm,"all")];
 LINEAR_REGRESSION = [nanmean(accuracy_test_linearreg);std(accuracy_test_linearreg);mean(time_linear_reg,"all")];
 NNET   = [nanmean(accuracy_test_nnet);std(accuracy_test_nnet);mean(time_net,"all")];
 SVM_RBF = [nanmean(accuracy_test_svmrbf);std(accuracy_test_svmrbf);mean(time_svm_rbf,"all")];
 NAIVE_BAYES= [nanmean(accuracy_test_nb);std(accuracy_test_nb);mean(time_nb,"all")];

 T = table(CNMF_SPA,CNMF_OPT,SPA_EM,CTD,MLE_AD,MLE_EM,SVM,LINEAR_REGRESSION,NNET,SVM_RBF,NAIVE_BAYES,...
    'RowNames',rows);

writetable(T,'table_car.csv','WriteRowNames',true)
save data_car