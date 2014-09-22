function CalcCorrespondence(w);

load('MeaMapPlot.mat','MeaMap');
% load('CorrDistributions2.mat', 'MaxCosSim')
load('CorrDistribution_Cult6_temp2.mat', 'LagAtMax', 'MaxCosSim')
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat', 'DataSet')
% MaxCosSim = MaxCosSim{w};
ic= DataSet{w}.ic;
a2a=MaxCosSim(size(ic,2)+1:end,size(ic,2)+1:end);
n2n=MaxCosSim(1:size(ic,2),1:size(ic,2));
a2n=MaxCosSim(1:size(ic,2),size(ic,2)+1:end);
ch = DataSet{w}.channel;
%---Some Top Connections---%
[~,~,~,a1,a_e2,~,ch,~]=ParseConnectionValues(w,1000);
astros=[];
for i=1:numel(unique(a_e2))
    [test(i),F(i),~]=mode(a_e2);
    astros = [astros,a1(ismember(a_e2,test(i)))];
    a_e2(a_e2==test(i))=[];
end

test(F<(mean(F)+std(F)))=[];
astros(F<(mean(F)+std(F)))=[];

%---ViewField Set of Electrodes---%
[~,~,channels] = FindElecsinVF(ch,MeaMap);
corrs = arrayfun(@(x) max(MaxCosSim(ismember(ic(1,:),channels),x)), ic(1,:));

%---Page Rank N2N---%
n2n(logical(isnan(n2n)))=0;
prn = pagerank(n2n);
prn = tiedrank(prn);
prn=prn';

%---Page Rank A2A---%
a2a(logical(isnan(a2a)))=0;
pra = pagerank(a2a);
pra = tiedrank(pra);
pra=pra';

%---Page Rank A2N---%
mat = zeros(size(a2n,1)+size(a2n,2),size(a2n,1)+size(a2n,2));
mat(1:size(ic,2),size(ic,2)+1:end)=MaxCosSim(1:size(ic,2),size(ic,2)+1:end);
mat=mat+mat';
mat(isnan(mat))=0;
pran = pagerank(mat);
prna = tiedrank(pran(1:size(ic,2)));
prna=prna';

pran = tiedrank(pran(size(ic,2)+1:end)); %ranks of Astrocytes in A2N
pran=pran';
%------%
ingroup = find(ismember(ic(1,:),test));
outgroup = setdiff(1:numel(ic(1,:)),find(ismember(ic(1,:),test)));

ingroupa = find(ismember(unique(a1),unique(astros)));
outgroupa = setdiff(1:numel(unique(a1)),ingroupa);

subplot(1,4,1)
title('Correlations');
notBoxPlot([padarray(corrs(ingroup),[0,numel(outgroup)-numel(ingroup)],nan,'post');... %corr of vf elecs to main pairs
            corrs(outgroup)]'); % corr of vf elecs to rest
set(gca,'ylim',[0,1],'XTickLabel', {'VF Elecs vs. main pairs', 'VF Elecs vs. rest'});
XYrotalabel(45,0);

subplot(1,4,2:4)
title('Centrality');  
notBoxPlot([padarray(prna(channels),[0,numel(outgroup)-numel(channels)],nan,'post');... %page rank of vf elecs in a2n
            prna(outgroup);... %page rank of rest in a2n
            padarray(prna(ingroup),[0,numel(outgroup)-numel(ingroup)],nan,'post');...   %page rank of main elecs in a2n
            padarray(prn(ingroup),[0,numel(outgroup)-numel(ingroup)],nan,'post');...    %page rank of main elecs in n2n
            padarray(prn(channels),[0,numel(outgroup)-numel(channels)],nan,'post');...  %page rank of vf elecs in n2n
            prn(outgroup);...                                                           %page rank of rest of elecs in n2n
            padarray(pra(ingroupa),[0,numel(outgroup)-numel(ingroupa)],nan,'post');...   %page rank of astro's connected to main elecs
            padarray(pra(outgroupa),[0,numel(outgroup)-numel(outgroupa)],nan,'post')]');  %page rank of astro's connected to other elecs                                                   
ylabel('Page Rank (a.u.)');
set(gca,'XTickLabel', {'A2N VF Elecs','A2N Rest of Elecs', 'A2N Main Elecs', 'N2N Main Elecs','N2N VF Elecs','N2N Rest of Elecs','A2A Main Elecs','A2A Rest of Elecs'});
XYrotalabel(45,0);
end