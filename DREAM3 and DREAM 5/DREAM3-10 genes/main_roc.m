%compute the AUC of epc_***
clear

%% choose data
%Input:
%(1)data: ngene*nsample; (2)G_gold: ture adjacency matrix
%(3)corr_G: correlation matrix; (4)mi_G: mutual information matrix

load('10gene.mat');

n_gene = size(data, 1);

%% choose lambda
lambda = 0.001:0.001:0.130; %lambda for CMI, MPMI, PMI and NPA
%% auc
L = length(lambda);
TPR = zeros(L, 1); FPR = zeros(L, 1);
TP = zeros(L, 1); FP = zeros(L, 1); 
FN = zeros(L, 1); TN = zeros(L, 1);
for i = 1:length(lambda)
    G = epc_cmi_vv1(mi_G, corr_G, lambda(i)); %CMI
%   G = epc_nonlinear_vv1(mi_G, corr_G, lambda(i)); %NPA
%    G = epc_mpmi_vv1(mi_G, corr_G, lambda(i)); % MPMI
%   G = epc_pmi_vv1(mi_G, corr_G, lambda(i)); % PMI

        
    G = G + tril(2*ones(size(G))); %set diag and lowtriange elements to 2
    [r1, c1] = find(G ==1); %get the position of positive edges (=1) in G
    temp_v1 = G_gold(sub2ind(size(G), r1, c1)); %get the true value of positive edges
    [r2, c2] = find(G ==0); %get the position of negative edges (=0) in G
    temp_v2 = G_gold(sub2ind(size(G), r2, c2)); %get the true value of negative edges
    TP(i) = sum(temp_v1); FP(i) = length(temp_v1) - TP(i);
    FN(i) = sum(temp_v2); TN(i) = length(temp_v2) - FN(i);
    TPR(i) = TP(i)/(TP(i)+FN(i)); FPR(i) = FP(i)/(FP(i) + TN(i));
    sprintf('%d/%d finished!\n', i, L)
end
roc_curve = unique([0, 0; FPR, TPR; 1, 1], 'rows');
auc = trapz(roc_curve(:, 1), roc_curve(:, 2));
