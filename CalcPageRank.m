function [val,rank] = CalcPageRank(mat);

Centrality = pagerank(mat+mat');
Centrality = (Centrality - min(Centrality))./(max(Centrality)-min(Centrality));
[val,rank]=sort(Centrality,'descend');
rank(isnan(val))=nan;

end