function [G, stepflag]  = epc_nonlinear_vv1(mi_G, corr_G, lambda)
%EPC_NONLINEAR_VV1 Compute the network structure using extended path consistent
%algorithm with nonlinear index NPA as the indicator for direct connection.
%Gaussian assumption is used.

%Input:
%mi_G: ngene*ngene mutual information matrix
%corr_G: the correlation matrix
%lambda: Threshold for 0\1 edges' distinguishment
%
%Output:
%G: ngene*ngene adject matrix, only uptriangle part has nonzero elements.
%stepflag: the iteration steps number.

G=(mi_G > lambda);

[edgerow, edgecol] = find(G); edgecount = numel(edgerow);

stop = (edgecount==0); Nstep = 10; stepflag = 0;
while (~stop)
    stepflag = stepflag + 1;
    %logical vectors of adjececency nodes
    indx = cellfun( @(r, c) (G(r, :)|G(:, r)'|G(c, :)|G(:, c)'), num2cell(edgerow), num2cell(edgecol), 'UniformOutput', 0);
    %exclude the nodes themselves from the adjecency nodes vectors
    indx = cellfun(@(x, r, c) indxexclude(x, r, c), indx, num2cell(edgerow), num2cell(edgecol), 'UniformOutput', 0);
    %find the non-all-zero logical vector
    indx2 = cellfun(@(x) any(x), indx);
    %find the indexes of adjecency nodes
    indx22 = cellfun(@(x) find(x), indx(indx2), 'UniformOutput', 0);
    tempG = cellfun(@(r, c, v) (mscmi(corr_G([r, c, v], [r, c, v]))>lambda), num2cell(edgerow(indx2)), num2cell(edgecol(indx2)), indx22);
    
    indx3 = sub2ind(size(G), edgerow(indx2), edgecol(indx2)); % the edges need to change
    stop = (norm(G(indx3) - tempG, Inf)==0) | (stepflag > Nstep);
    G(indx3) = tempG;
    [edgerow, edgecol] = find(G);
end

end

%% compute the indx position excluding x,y themselves
function out = indxexclude(x, r, c)
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
cmiv = cmi0*Covm1^2/Covm3/Covm2; %CMI(x,y|z)*|C(z)|^2/|C(x,z)|/|C(y,z)|

% cmiv=abs(cmiv);
if  cmiv==inf 
    cmiv=0;
end
end
