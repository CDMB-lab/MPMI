%compute the CMI,NPA,MPMI,PMI with Gaussian assumption

function out = scalecmi(epsilon)
L = 5;
%% choose strong connection matrix
%B = [15,8,1;8,15,1;1,1,15] + 1/epsilon*[1,0,-1;0,0,0;-1,0,1]; %one strong connection
B = [15,8,1;8,15,1;1,1,15] + 1/epsilon*[1,0,-1;0,1,-1;-1,-1,2]; %two strong connections

%% computing indexes
%generate data
A = inv(B);
data = mvnrnd([0, 0, 0], A, 10^6);
out = zeros(L,1);
D = corr(data);

out_patrtial = abs(partialcorr(data(:, 1), data(:, 2), data(:, 3)));
out(1) = -0.5*log(1-out_patrtial^2); %CMI for Gaussian or near-Gaussian
% out(2) = exp(2*mi_new(data(:, 1), data(:, 3)))*out(2)*exp(2*mi_new(data(:, 2), data(:, 3))); %NPA
out(3) = out(1)/(1-corr(data(:, 1), data(:, 3))^2)/(1-corr(data(:,2), data(:, 3))^2);%NPA for Gaussian
out(4) = mpmi(D);%MPMI for Gaussian
out(5) = pmi(D);%PMI for Gaussian

end
function mpmiv=mpmi(D)
    n = size(D,1);
    Covm = D;
    Cov1 = Covm(1,1);%c(x)
    Cov2 = Covm(2,2);%c(y)
    Covm1 = Covm(3:end,3:end);%c(z)
    Covm2 =Covm([1, 3:end], [1, 3:end]);%c(x,z)
    Covm3 = Covm([2, 3:end], [2, 3:end]);%c(y,z)

    InvCov1 = 1/Cov1;
    InvCov2 = 1/Cov2;
    InvCovm = inv(Covm);
    InvCovm1 = inv(Covm1);
    InvCovm2 = inv(Covm2);
    InvCovm3 = inv(Covm3);

    C11 = -InvCovm(1,2)*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*InvCovm(1,2)+InvCovm(1,1);
    C12 = 0;
    C13 = -InvCovm(1,2)*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:end)-InvCovm3(1,2:end))+InvCovm(1,3:end);
    C23 = -InvCovm(1,2)*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:end)-InvCovm2(1,2:end))+InvCovm(2,3:end);
    C22 = -InvCovm(1,2)*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*InvCovm(1,2)+InvCovm(2,2);
    C33 = -(InvCovm(2,3:end)-InvCovm3(1,2:end))'*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:end)-InvCovm3(1,2:end))-(InvCovm(1,3:end)-InvCovm2(1,2:end))'*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:end)-InvCovm2(1,2:end))+(InvCovm(3:end,3:end)-InvCovm3(2:end,2:end))+(InvCovm(3:end,3:end)-InvCovm2(2:end,2:end))+InvCovm1;
    InvC = [[C11,C12,C13];[C12,C22,C23];[[C13',C23'],C33]];
    % C = inv(InvC);

    C0 = (det(Covm)*det(Covm1)/(det(Covm2)*det(Covm3)))*Cov1*Cov2*(InvCovm(2,2)-InvCovm3(1,1)+InvCov2)*(InvCovm(1,1)-InvCovm2(1,1)+InvCov1);
    pmiv = 0.5 * (trace(InvC*Covm)+log(C0)-n);

    cmi0 = 0.5*log(det(Covm2)*det(Covm3)/det(D)/det(Covm1));
    cmi1 = (det(Covm1)^2-det(Covm3)*det(Covm2))/det(Covm3)/det(Covm2); 

    mpmiv = pmiv + cmi0*cmi1;
end
function pmiv=pmi(D)
    n = size(D,1);
    Covm = D;
    Cov1 = Covm(1,1);%c(x)
    Cov2 = Covm(2,2);%c(y)
    Covm1 = Covm(3:end,3:end);%c(z)
    Covm2 =Covm([1, 3:end], [1, 3:end]);%c(x,z)
    Covm3 = Covm([2, 3:end], [2, 3:end]);%c(y,z)

    InvCov1 = 1/Cov1;
    InvCov2 = 1/Cov2;
    InvCovm = inv(Covm);
    InvCovm1 = inv(Covm1);
    InvCovm2 = inv(Covm2);
    InvCovm3 = inv(Covm3);

    C11 = -InvCovm(1,2)*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*InvCovm(1,2)+InvCovm(1,1);
    C12 = 0;
    C13 = -InvCovm(1,2)*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:end)-InvCovm3(1,2:end))+InvCovm(1,3:end);
    C23 = -InvCovm(1,2)*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:end)-InvCovm2(1,2:end))+InvCovm(2,3:end);
    C22 = -InvCovm(1,2)*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*InvCovm(1,2)+InvCovm(2,2);
    C33 = -(InvCovm(2,3:end)-InvCovm3(1,2:end))'*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:end)-InvCovm3(1,2:end))-(InvCovm(1,3:end)-InvCovm2(1,2:end))'*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:end)-InvCovm2(1,2:end))+(InvCovm(3:end,3:end)-InvCovm3(2:end,2:end))+(InvCovm(3:end,3:end)-InvCovm2(2:end,2:end))+InvCovm1;
    InvC = [[C11,C12,C13];[C12,C22,C23];[[C13',C23'],C33]];
    % C = inv(InvC);

    C0 = (det(Covm)*det(Covm1)/(det(Covm2)*det(Covm3)))*Cov1*Cov2*(InvCovm(2,2)-InvCovm3(1,1)+InvCov2)*(InvCovm(1,1)-InvCovm2(1,1)+InvCov1);
    pmiv = 0.5 * (trace(InvC*Covm)+log(C0)-n);
end
