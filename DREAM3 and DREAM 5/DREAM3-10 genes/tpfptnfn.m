function out = tpfptnfn(G, G_gold)
out = zeros(9, 1);

G = G + tril(2*ones(size(G))); %set diag and lowtriange elements to 2
[r1, c1] = find(G ==1); %get the position of positive edges (=1) in G
temp_v1 = G_gold(sub2ind(size(G), r1, c1)); %get the true value of positive edges
[r2, c2] = find(G ==0); %get the position of negative edges (=0) in G
temp_v2 = G_gold(sub2ind(size(G), r2, c2)); %get the true value of negative edges
out(1) = sum(temp_v1); %TP
out(2) = length(temp_v1) - out(1); %FP
out(3) = sum(temp_v2); %FN
out(4) = length(temp_v2) - out(3); %TN
out(5) = out(1)/(out(1) + out(3)); %TPR
out(6) = out(2)/(out(2) + out(4)); %FPR
out(7) = out(1)/(out(1) + out(2)); %Precision
out(8) = out(4)/(out(3) + out(4)); %NPV
out(9) = (out(1)+out(4))/(out(1)+out(2)+out(3)+out(4)); %ACC

end