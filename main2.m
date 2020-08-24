clear all

predictionMethod = 'ppxa'%

%----add dependencies to path----
addpath(genpath('helper_functions'));

%----read data---- % source=HNRD, another: https://github.com/cuizhensdws/drug-disease-datasets
data_dir='data/'
datasets={'Fdatasets', 'Cdatasets'}

%----define parameters----
n = 10;% 'n' in "n-fold experiment"
%global f_roc f_pr

for ds=[ 1 2 ]
%f_roc = figure;
%f_pr = figure;
   
dataname=datasets{ds};

%read data
load([data_dir dataname '/DiDrA.txt']);  load([data_dir dataname '/DiseaseSim.txt']);  load([data_dir dataname '/DrugSim.txt']);
Y=DiDrA;  Sd=DiseaseSim; St=DrugSim;

len=length(find(Y==1)); %or sum(sum(Y));
true_association_indices=find(Y==1);

rng(0);
rand_ind = true_association_indices(randperm(len));%randperm(len);

global pp mu1 mu2 lamda
getParameters(predictionMethod,ds)
        
    % loop over the n folds
    AUCs  = zeros(1,n);  AUPRs = zeros(1,n); 
    XsROC =[]; YsROC =[]; XsPR =[]; YsPR =[];
    pr=[]; re=[]; acc=[]; F=[];
    for i=1:n

       test_ind = rand_ind((floor((i-1)*len/n)+1:floor(i*len/n))');
       test_ind = test_ind(:);
       
       
       y2 = Y;
       y2(test_ind) = 0;
       fprintf('*');

       tic
       y3 = eval([ predictionMethod  '(y2,Sd,St,test_ind)']);
       
       time_taken=toc
       
       test_ind2=find(y2==0);
       test_ind2=test_ind2(:);
       
       
       [AUCs(i),XcROC,YcROC]  = calculate_auc (y3(test_ind2),Y(test_ind2), 0);    
       [AUPRs(i),XcPR,YcPR, T] = calculate_aupr(y3(test_ind2),Y(test_ind2), 0);
     
    end
    
    auc= mean(AUCs) %round( mean(AUCs), 2)
    aupr=mean(AUPRs)%round( mean(AUPRs), 2)
   
end

