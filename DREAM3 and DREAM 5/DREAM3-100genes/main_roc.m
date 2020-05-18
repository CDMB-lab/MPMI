%compute the ROC of epc_***
clear

%% choose data
%With at least two matrices in the input:
%(1)data: ngene*nsample; (2)G_gold: ture adjacency matrix
%And two other matrices useful:
%(3)corr_G: correlation matrix; (4)mi_G:
%mutual information matrix


load('100gene.mat');


n_gene = size(data, 1);
% N = n_gene*(n_gene-1)/2;

%FOR 100 gene
% lambda = 10^(-5)*1.22.^(0:1:60); %lambda for CMI, MCMI, p12, nonlinear
% lambda = 0.001+0.001*1.1.^(0:1:60); %lambda for PMI
lambda = 0.001:0.001:0.130;


%FOR dream5 data 4511 gene / insilico 1643 gene
% lambda = 10^(-5)*1.22.^(10:1:50);
% lambda = fliplr(lambda);

%For 10^3gene with p=0.01 edges
% lambda = 0.01:0.01:1;
% lambda = 0.1:0.1:1;

%% compute roc choose algorithm

% %compute the correlation matrix if not included in the input data
% corr_G = corrcoef(data');
% %compute MI matrix if not included in the input data
% % mi_G = tril(ones(n_gene), -1)';
% % [r, c] = find(mi_G);
% % indx = sub2ind(size(mi_G), r, c);
% % tempG = arrayfun(@(i, j) -0.5*log(1-corr_G(i, j).^2), r, c);
% % mi_G(indx) = tempG;
% mi_G =  -0.5*log(1-corr_G.^2);
% mi_G = tril(mi_G, -1)';

L = length(lambda);
TPR = zeros(L, 1); FPR = zeros(L, 1);
TP = zeros(L, 1); FP = zeros(L, 1); 
FN = zeros(L, 1); TN = zeros(L, 1);
for i = 1:length(lambda)
%     G = epc_cmi_vv1(mi_G, corr_G, lambda(i)); %CMI
     G = epc_mpmi_vv1(mi_G, corr_G, lambda(i)); % MPMI
%      G = epc_pmi_vv1(mi_G, corr_G, lambda(i)); % PMI
%     G = epc_p12_vv1(mi_G, corr_G, lambda(i)); %LPA
%    G = epc_nonlinear_vv1(mi_G, corr_G, lambda(i)); %NPA

        
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

% figure; plot(roc_curve(:, 1), roc_curve(:, 2));
% xlabel('FPR'); ylabel('TPR')
% title('ROC curve')

% [x, idx] = unique(roc_curve(:, 1)); y = roc_curve(idx, 2);
% auc = sum( (y(1:end-1)+y(2:end))/2.*diff(x) );
auc = trapz(roc_curve(:, 1), roc_curve(:, 2));

%% save different data
% save('100_cmi_roc.mat', 'roc_curve', 'auc');
% save('100_mcmi_roc.mat', 'roc_curve', 'auc');
% save('100_p12_roc.mat', 'roc_curve', 'auc');
% save('100_nonlinear_roc.mat', 'roc_curve', 'auc');

% save('d5_cmi_roc.mat', 'roc_curve', 'auc');
% save('d5_mcmi_roc.mat', 'roc_curve', 'auc');
% save('d5_p12_roc.mat', 'roc_curve', 'auc');
% save('d5_nonlinear_roc.mat', 'roc_curve', 'auc');

% save('si_cmi_roc.mat', 'roc_curve', 'auc', 'TP', 'FP', 'FN', 'TN', 'TPR', 'FPR');
% save('si_mcmi_roc.mat', 'roc_curve', 'auc', 'TP', 'FP', 'FN', 'TN', 'TPR', 'FPR');
% save('si_p12_roc.mat', 'roc_curve', 'auc', 'TP', 'FP', 'FN', 'TN', 'TPR', 'FPR');
% save('si_non_roc.mat', 'roc_curve', 'auc', 'TP', 'FP', 'FN', 'TN', 'TPR', 'FPR');

