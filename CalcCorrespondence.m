function CalcCorrespondence(w);

load('MeaMapPlot.mat','MeaMap');
load('CorrDistributions2.mat', 'MaxCosSim')
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat', 'DataSet')
MaxCosSim = MaxCosSim{w};
ic= DataSet{w}.ic;
a2a=MaxCosSim(size(ic,2)+1:end,size(ic,2)+1:end);
n2n=MaxCosSim(1:size(ic,2),1:size(ic,2));
a2n=MaxCosSim(1:size(ic,2),size(ic,2)+1:end);
ch = DataSet{w}.channel;
%---Some Top Connections---%
[~,~,~,~,a_e2,aval,ch,~]=ParseConnectionValues(w,1000);
for i=1:numel(unique(a_e2))
    [test(i),~,~]=mode(a_e2);
    a_e2(a_e2==test(i))=[];
end

%---ViewField Set of Electrodes---%
[~,~,channels] = FindElecsinVF(ch,MeaMap);
corrs = arrayfun(@(x) max(MaxCosSim(ismember(ic(1,:),channels),x)), ic(1,:));

%---Page Rank---%
mat = zeros(size(a2n,1)+size(a2n,2),size(a2n,1)+size(a2n,2));
mat(1:size(ic,2),size(ic,2)+1:end)=MaxCosSim(1:size(ic,2),size(ic,2)+1:end);
mat=mat+mat';
mat(isnan(mat))=0;
pr = pagerank(mat);
pr = tiedrank(pr(1:size(ic,2)));
pr=pr';

ingroup = find(ismember(ic(1,:),test));
outgroup = setdiff(1:numel(ic(1,:)),find(ismember(ic(1,:),test)));
subplot(1,2,1)
notBoxPlot([padarray(corrs(ingroup),[0,numel(outgroup)-numel(ingroup)],nan,'post');corrs(outgroup)]');
set(gca,'ylim',[0,1],'XTickLabel', {'VF Elecs vs. main pairs', 'VF Elecs vs. rest'});
ylabel('Correlation');
subplot(1,2,2)
notBoxPlot([padarray(pr(ingroup),[0,numel(outgroup)-numel(ingroup)],nan,'post');pr(outgroup)]');
set(gca,'XTickLabel', {'A2N: Main Connections', 'A2N: Rest of Connections'});
ylabel('Page Rank');
[h1,p1] = ttest2(corrs(ingroup),corrs(outgroup))
[p1,h1] = ranksum(pr(ingroup),pr(outgroup))


end