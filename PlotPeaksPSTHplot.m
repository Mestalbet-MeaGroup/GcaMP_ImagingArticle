function PlotPeaksPSTHplot(CutoutNear,ranksC,fs,minx,maxx);
mat = CutoutNear(ranksC==3|ranksC==4,:)';
test=pdist(corrcoef(zscore(mat,[],2)),'spearman');
tree=linkage(test,'ward');
ClusterIDs  = cluster(tree,'MaxClust',10,'Criterion','distance');
hold on
% color = cbrewer('qual','Set2',max(ClusterIDs));
color = bone(max(ClusterIDs)+5);
for i=1:max(ClusterIDs)
    %     subplot(floor(sqrt(max(ClusterIDs))),ceil(sqrt(max(ClusterIDs))),i);
    plotmat = mean(mat(:,ClusterIDs==i),2);
    [~,pkloc]=max(plotmat);
    plot(([1:size(plotmat,1)]-pkloc)./fs,plotmat,'.-','Color',color(i,:));   
end
axis tight;
xlim([minx,maxx]);
axis off;
end