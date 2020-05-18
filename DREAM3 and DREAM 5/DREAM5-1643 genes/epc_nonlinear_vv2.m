function [G, stepflag]  = epc_nonlinear_vv2(mi_G, corr_G, lambda)
%EPC_NONLINEAR_VV2 Compute the network structure using extended path consistent
%algorithm with NPA as the indicator for direct connection. 
%Gaussian assumption is used.

%Input:
%mi_G: ngene*ngene mutual information matrix
%corr_G: the correlation matrix
%lambda: Threshold for 0\1 edges' distinguishment
%
%Output:
%G: ngene*ngene adject matrix, only uptriangle part has nonzero
%elements.
%stepflag: the iteration steps number.

G=(mi_G > lambda);

[edgerow, edgecol] = find(G); edgecount = numel(edgerow);
mi_G=mi_G+(mi_G)';
stop = (edgecount==0); Nstep = 5; stepflag = 0;
while (~stop)
    stepflag = stepflag + 1;
    %logical vectors of adjececency nodes
    a1 = cellfun( @(r) (G(r, :)|G(:, r)'), num2cell(edgerow), 'UniformOutput', 0);
    a2 = cellfun( @(r) (G(r, :)|G(:, r)'), num2cell(edgecol), 'UniformOutput', 0);
    a1 = cellfun(@(x, c) indxexclude(x, c), a1, num2cell(edgecol), 'UniformOutput', 0);
    a2 = cellfun(@(x, c) indxexclude(x,c), a2, num2cell(edgerow), 'UniformOutput', 0);
    b1 = ones(size(a1,1),1);
    b2 = ones(size(a2,1),1);
    e1= cellfun(@(x) find(x),a1(logical(b1)), 'UniformOutput', 0);
    e2= cellfun(@(x) find(x),a2(logical(b2)), 'UniformOutput', 0);
    temp1 = cellfun(@(r,v) mi_G(r,v), num2cell(edgerow),e1,'UniformOutput',false);
    temp2 = cellfun(@(r,v) mi_G(r,v), num2cell(edgecol),e2,'UniformOutput',false);
    [~,ix_a1]=cellfun(@(r) sort(r,'descend'), temp1, 'UniformOutput',false);
    [~,ix_a2]=cellfun(@(r) sort(r,'descend'), temp2, 'UniformOutput',false);
    s1=cellfun(@(r,c) r(1:num(c)), ix_a1,e1,'UniformOutput',false);
    s2=cellfun(@(r,c) r(1:num(c)), ix_a2,e2, 'UniformOutput',false);
    f1=cellfun(@(r,c) r(c), e1,s1, 'UniformOutput',false);
    f2=cellfun(@(r,c) r(c), e2,s2, 'UniformOutput',false);
    indx22 = cellfun(@(r,c) union(r,c), f1,f2, 'UniformOutput',false);
    tempG = cellfun(@(r, c, v) (npa(corr_G([r, c, v], [r, c, v]))>lambda), num2cell(edgerow), num2cell(edgecol), indx22);
    
    indx3 = sub2ind(size(G), edgerow, edgecol); % the edges need to change
    stop = (norm(G(indx3) - tempG, Inf)==0) | (stepflag > Nstep);
    G(indx3) = tempG;
    [edgerow, edgecol] = find(G);
end

end

%% compute the indx position excluding x,y themselves
function out = indxexclude(x,c)
 x(c) = 0;
out = x;
end
%% choose the most correlated five variables
function out = num(x)
a = numel(x);
if a>5
    out = 5;
else
    out = a;
end
end
%% compute CMI of x and y
function npav=npa(corr_G)
%the first two coordinates are x and y, others are adjecency nodes
Covm = det(corr_G); %C(x,y,z)
Covm1 = det(corr_G(3:end, 3:end)); %C(z)
Covm2 = det(corr_G([1, 3:end], [1, 3:end])); %C(x,z)
Covm3 = det(corr_G([2, 3:end], [2, 3:end])); %C(y,z)

cmi0 = 0.5*log(Covm2*Covm3/Covm/Covm1);

npav = cmi0*Covm1^2/Covm3/Covm2; %CMI(x,y|z)*|C(z)|^2/|C(x,z)|/|C(y,z)|

% cmiv=abs(cmiv);
if  npav==inf 
    npav=0;
end
end
