function [bg, indegree,outdegree,dist,reachable,emptynodes,noderank,edgerank]=CalcAstroGraph(a2a,lagmat,cutoff)
%% Prepare graph adjacency matrix
% cutoff = triu(a2a,1)+tril(nan(size(a2a)));
% cutoff = nanmean(cutoff(:))+nanstd(cutoff(:));
% cutoff = 0;
conmatWithWeights = (triu(a2a,1)+triu(a2a,1)');
lagmat = (triu(lagmat,1)+triu(lagmat,1)');
conmatWithWeights(conmatWithWeights<cutoff)=0;
lagmat(conmatWithWeights<cutoff)=0;

emptynodes = find(sum([sum(conmatWithWeights,1)==0;[sum(conmatWithWeights,2)==0]'],1)==2);
conmatWithWeights(emptynodes,:)=[];
conmatWithWeights(:,emptynodes)=[];
lagmat(emptynodes,:)=[];
lagmat(:,emptynodes)=[];
lagmat=lagmat+lagmat';
conmatWithWeights(triu(lagmat)>0)=0;
conmatWithWeights(tril(lagmat)<0)=0;
%% Calculate node centrality 

bg = biograph(sparse(conmatWithWeights));

Centrality = pagerank(conmatWithWeights>0);
[~,~,noderank] = unique(Centrality);
val = conmatWithWeights(conmatWithWeights>0);
[~,~,edgerank] = unique(round(val*1000)/1000,'stable');
edgecolor = parula(max(edgerank));

set(bg.Nodes,'Shape','circle');
for i=1:size(bg.edges,1)
    set(bg.edges(i),'LineWidth',val(i)*1.5);
    set(bg.edges(i),'LineColor',edgecolor(edgerank(i),:));
end

nodecolor = parula(max(noderank));
for i=1:size(bg.Nodes,1)
    set(bg.Nodes(i),'Color',nodecolor(noderank(i),:));
    set(bg.Nodes(i),'Size',[noderank(i),noderank(i)]);
    set(bg.Nodes(i),'ID',num2str(i));
    set(bg.Nodes(i),'FontSize',24);
end
set(bg,'ShowArrows','on','ShowWeights','off','NodeAutoSize','off')
% view(bg)
%% Calculate network statistics

[dist] = allshortestpaths(bg,'Directed',true);
indegree = sum(conmatWithWeights);
outdegree = sum(conmatWithWeights,2)';
[reachable, ~] = conncomp(bg,'Directed',true);
end