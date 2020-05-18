%compute the rightlambda of epc_***
clear

%% choose data
%Input:
%(1)data: ngene*nsample; (2)G_gold: ture adjacency matrix
%(3)corr_G: correlation matrix; (4)mi_G:
%mutual information matrix

load('10gene.mat');

n_gene = size(data, 1);

%% choose lambda
%for 10 and 100 genes
lambda = 0.001:0.001:0.130;%lambda for CMI, MCMI, LPA, NPA
%for 1643 genes
% lambda = 10^(-4)*1.22.^(10:1:50); %lambda for CMI, MCMI, LPA, NPA
%% compute roc choose algorithm

L = length(lambda);

L_index = 9; %number of indexes in roc
roc_index = zeros(L, L_index);
for i = 1:length(lambda)
%   G = epc_cmi_vv1(mi_G, corr_G, lambda(i)); %CMI
    G = epc_mpmi_vv1(mi_G, corr_G, lambda(i)); %MPMI
%   G = epc_pmi_vv1(mi_G, corr_G, lambda(i)); %PMI
%   G = epc_nonlinear_vv1(mi_G, corr_G, lambda(i)); %NPA

    roc_index(i, :) = tpfptnfn(G, G_gold);
    
    sprintf('%d/%d finished!\n', i, L)
end

%choose best lambda
FPR = roc_index(:, 6); TPR = roc_index(:, 5); ACC = roc_index(:, 9);
[~, lambdainx] = max(ACC);
% [~, lambdainx] = max(roc_index(:, 7)+roc_index(:, 8));%7 Precision=TP/(TP+FP),8 Negative predictive value=TN/(FN+TN), 
rightlambda = lambda(lambdainx);
