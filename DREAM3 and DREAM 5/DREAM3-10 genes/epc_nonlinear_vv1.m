function [G, stepflag]  = epc_nonlinear_vv1(mi_G, corr_G, lambda)
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

G=(mi_G > lambda);%大于阈值的为1，小于的为0,此时G为逻辑矩阵

[edgerow, edgecol] = find(G); edgecount = numel(edgerow);%find函数查询矩阵G中非零元素的位置

stop = (edgecount==0); Nstep = 10; stepflag = 0;
while (~stop)
    stepflag = stepflag + 1;
    %logical vectors of adjececency nodes邻接节点的逻辑向量,每个节点都对应着一个逻辑向量
    indx = cellfun( @(r, c) (G(r, :)|G(:, r)'|G(c, :)|G(:, c)'), num2cell(edgerow), num2cell(edgecol), 'UniformOutput', 0);
    %exclude the nodes themselves from the adjecency nodes
    %vectors从邻接节点向量中排除节点,不清楚是如何排除的
    indx = cellfun(@(x, r, c) indxexclude(x, r, c), indx, num2cell(edgerow), num2cell(edgecol), 'UniformOutput', 0);
    %find the non-all-zero logical vector找到非全零的逻辑向量
    indx2 = cellfun(@(x) any(x), indx);
    %find the indexes of adjecency nodes查找邻接节点的索引
    indx22 = cellfun(@(x) find(x), indx(indx2), 'UniformOutput', 0);
    tempG = cellfun(@(r, c, v) (mscmi(corr_G([r, c, v], [r, c, v]))>lambda), num2cell(edgerow(indx2)), num2cell(edgecol(indx2)), indx22);
    
    indx3 = sub2ind(size(G), edgerow(indx2), edgecol(indx2)); % the edges need to change
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

%% compute CMI of x and y
function cmiv=mscmi(corr_G)
%the first two coordinates are x and y, others are adjecency nodes
Covm = det(corr_G); %C(z,y,z)
Covm1 = det(corr_G(3:end, 3:end)); %C(z)
Covm2 = det(corr_G([1, 3:end], [1, 3:end])); %C(y,z)
Covm3 = det(corr_G([2, 3:end], [2, 3:end])); %C(x,z)

cmi0 = 0.5*log(Covm2*Covm3/Covm/Covm1);

% cmiv = cmi0/Covm3/Covm2; %CMI(x,y|z)/|C(x,z)|/|C(y,z)|
%As |C(x)| = |C(y)| =1
cmiv = cmi0*Covm1^2/Covm3/Covm2; %CMI(x,y|z)*|C(z)|^2/|C(x,z)|/|C(y,z)|

% cmiv=abs(cmiv);
if  cmiv==inf 
    cmiv=0;
end
end
