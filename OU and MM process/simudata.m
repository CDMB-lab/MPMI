%simulate data by OU and MM process model

clear

load('test3.mat');
Ngene = size(G_gold, 1);

deltat = 10^-5; N_iter0 = 10^5; sigma =4^2;%for mm process %sigma =0.4;%for ou process
N_data = 10^6;

A = G_gold;
epsilon = 4^(-5);

LAMbda = zeros(size(G_gold)); LAMbda(sub2ind(size(G_gold), r(1:2), c(1:2))) = -1;
LAMbda = LAMbda + LAMbda';
LAMbda = LAMbda - diag(sum(LAMbda, 1));

B = 2*A +diag(sum(2*A)+10)+ LAMbda/epsilon;

%  MM model add the next two line 
    B0 = diag(B);
    B = tril(B, -1); B = B + B';

    
    x = zeros(Ngene, 1);
    for j = 1:N_iter0
        %MM model
        x = x - deltat*(B0.*x + B*(x./(1+10^(-4)*x))) + sqrt(deltat*2*sigma)*randn(Ngene, 1);
        %OU model
%       x = x - deltat*(B*x) + sqrt(deltat*2*sigma)*randn(Ngene, 1);
    end
    data = zeros(Ngene, N_data);
    data(:, 1) = x;
    for j = 2:N_data
        for k = 1:10
            %MM model
            x = x - deltat*(B0.*x + B*(x./(1+10^(-4)*x))) + sqrt(deltat*2*sigma)*randn(Ngene, 1);
            %OU model
%           x = x - deltat*(B*x) + sqrt(deltat*2*sigma)*randn(Ngene, 1);
        end
        data(:, j) = x;
    end
    
    corr_G = corr(data');
    mi_G = -0.5*log(1-corr_G.^2);
    mi_G = tril(mi_G, -1)';
    
save('mmdata_acc.mat', 'corr_G', 'mi_G', 'G_gold', 'data');
% save('oudata_acc.mat', 'corr_G', 'mi_G', 'G_gold', 'data');
        
