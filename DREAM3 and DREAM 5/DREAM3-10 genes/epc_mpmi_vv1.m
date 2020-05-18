function [G, stepflag]  = epc_mpmi_vv1(mi_G, corr_G, lambda)
%EPC_NONLINEAR_VV1 Compute the network structure using extended path consistent
%algorithm with NPA = EXP(MI(x, z))*CMI(x, y|z)*EXP(MI(y, z)) as 
%the indicator for direct connection. Gaussian assumption is used.

%Input:
%mi_G: ngene*ngene mutual information matrix
%corr_G: the correlation matrix
%lambda: Threshold for 0\1 edges' distinguishment
%
%Output:
%G: ngene*ngene adject matrix, only uptriangle part has nonzero
%elements.%应该是指邻接矩阵，只有上三角部分有非零元素
%stepflag: the iteration steps number.迭代步骤数目

G=(mi_G > lambda);%大于阈值的为1，小于的为0

[edgerow, edgecol] = find(G); edgecount = numel(edgerow);%find函数查询矩阵G中非零元素的位置

stop = (edgecount==0); Nstep = 10; stepflag = 0;
while (~stop)
    stepflag = stepflag + 1;
    %logical vectors of adjececency nodes邻接节点的逻辑向量
    indx = cellfun( @(r, c) (G(r, :)|G(:, r)'|G(c, :)|G(:, c)'), num2cell(edgerow), num2cell(edgecol), 'UniformOutput', 0);
    %exclude the nodes themselves from the adjecency nodes vectors从邻接节点向量中排除节点
    indx = cellfun(@(x, r, c) indxexclude(x, r, c), indx, num2cell(edgerow), num2cell(edgecol), 'UniformOutput', 0);
    %find the non-all-zero logical vector找到非全零的逻辑向量
    indx2 = cellfun(@(x) any(x), indx);
    %find the indexes of adjecency nodes查找邻接节点的索引
    indx22 = cellfun(@(x) find(x), indx(indx2), 'UniformOutput', 0);
    tempG = cellfun(@(r, c, v) (mpmi(corr_G([r, c, v], [r, c, v]))>lambda), num2cell(edgerow(indx2)), num2cell(edgecol(indx2)), indx22);
    
    indx3 = sub2ind(size(G), edgerow(indx2), edgecol(indx2)); % the edges need to change找到G对应位置的数的序号
    stop = (norm(G(indx3) - tempG, Inf)==0) | (stepflag > Nstep);%norm(A,inf)即返回max(abs(A))
    G(indx3) = tempG;
    [edgerow, edgecol] = find(G);
end

end

%% compute the indx position excluding x,y themselves
function out = indxexclude(x, r, c)
% x(x ==r | x==c) = [];
% out = x;
x(r) = 0; x(c) = 0;
out = x;
end
%% compute MPMI of x and y
function mpmiv=mpmi(corr_G)
    n = size(corr_G,1);
    Covm = corr_G;
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
    C13 = -InvCovm(1,2)*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:end)-InvCovm3(1,2:end))+InvCovm(1,3:end);%2+n1即可以改为end,1+n1按理说也可以改为end
    C23 = -InvCovm(1,2)*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:end)-InvCovm2(1,2:end))+InvCovm(2,3:end);
    C22 = -InvCovm(1,2)*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*InvCovm(1,2)+InvCovm(2,2);
    % C33多加了一个InvCovm1即z的协方差矩阵的逆，同时里面有的地方用了转置来乘，公式中没有表现转置，可能因为这个原因多加了个z的逆
    C33 = -(InvCovm(2,3:end)-InvCovm3(1,2:end))'*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:end)-InvCovm3(1,2:end))-(InvCovm(1,3:end)-InvCovm2(1,2:end))'*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:end)-InvCovm2(1,2:end))+(InvCovm(3:end,3:end)-InvCovm3(2:end,2:end))+(InvCovm(3:end,3:end)-InvCovm2(2:end,2:end))+InvCovm1;
    InvC = [[C11,C12,C13];[C12,C22,C23];[[C13',C23'],C33]];
    % C = inv(InvC);

    C0 = (det(Covm)*det(Covm1)/(det(Covm2)*det(Covm3)))*Cov1*Cov2*(InvCovm(2,2)-InvCovm3(1,1)+InvCov2)*(InvCovm(1,1)-InvCovm2(1,1)+InvCov1);
    pmiv = 0.5 * (trace(InvC*Covm)+log(C0)-n);

    cmi0 = 0.5*log(det(Covm2)*det(Covm3)/det(corr_G)/det(Covm1));
    cmi1 = (det(Covm1)^2-det(Covm3)*det(Covm2))/det(Covm3)/det(Covm2); %实质计算的是公式后半部分

    mpmiv = pmiv + cmi0*cmi1;

    % cmiv=abs(cmiv);
    if  mpmiv==inf
        mpmiv=0;
    end
end
