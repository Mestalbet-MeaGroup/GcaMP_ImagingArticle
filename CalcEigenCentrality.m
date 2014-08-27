function [val,rank]=CalcEigenCentrality(mat,thresh)
%Accepts a symmetric matrix with and calculates the eigen centrality

test = triu(mat,1);
test(test==0)=nan;
% thresh = nanmean(test(:));
test(test>=thresh)=1;
test(test<thresh)=0;
test(isnan(test))=0;
test=test+test'+eye(size(test));
[V,D]=eig(test);
[max_eig,ind]=max(diag(D));
Centrality=V(:,ind);
[val,rank]=sort(Centrality,'descend');
end