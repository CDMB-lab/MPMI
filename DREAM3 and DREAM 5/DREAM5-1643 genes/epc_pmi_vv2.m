function [G, stepflag]  = epc_pmi_vv2(mi_G, corr_G, lambda)
%EPC_PMI_VV2 Compute the network structure using extended path consistent
%algorithm with PMI as the indicator for direct connection. Gaussian
%assumption is used.

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
    tempG = cellfun(@(r, c, v) (pmi(corr_G([r, c, v], [r, c, v]))>lambda), num2cell(edgerow), num2cell(edgecol), indx22);
    
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
%% compute MPMI of x and y
function pmiv=pmi(corr_G)
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
    C13 = -InvCovm(1,2)*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:end)-InvCovm3(1,2:end))+InvCovm(1,3:end);
    C23 = -InvCovm(1,2)*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:end)-InvCovm2(1,2:end))+InvCovm(2,3:end);
    C22 = -InvCovm(1,2)*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*InvCovm(1,2)+InvCovm(2,2);
    C33 = -(InvCovm(2,3:end)-InvCovm3(1,2:end))'*(1/(InvCovm(2,2)-InvCovm3(1,1)+InvCov2))*(InvCovm(2,3:end)-InvCovm3(1,2:end))-(InvCovm(1,3:end)-InvCovm2(1,2:end))'*(1/(InvCovm(1,1)-InvCovm2(1,1)+InvCov1))*(InvCovm(1,3:end)-InvCovm2(1,2:end))+(InvCovm(3:end,3:end)-InvCovm3(2:end,2:end))+(InvCovm(3:end,3:end)-InvCovm2(2:end,2:end))+InvCovm1;
    InvC = [[C11,C12,C13];[C12,C22,C23];[[C13',C23'],C33]];
    % C = inv(InvC);

    C0 = (det(Covm)*det(Covm1)/(det(Covm2)*det(Covm3)))*Cov1*Cov2*(InvCovm(2,2)-InvCovm3(1,1)+InvCov2)*(InvCovm(1,1)-InvCovm2(1,1)+InvCov1);
    pmiv = 0.5 * (trace(InvC*Covm)+log(C0)-n);

    % cmiv=abs(cmiv);
    if  pmiv==inf
        pmiv=0;
    end
end

