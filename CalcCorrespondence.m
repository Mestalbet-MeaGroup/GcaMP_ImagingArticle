function [a2a,a2n,n2n,MaxCosSim,mat,ic]=CalcCorrespondence(w);

load('MeaMapPlot.mat','MeaMap');
% load('CorrDistributions2.mat', 'MaxCosSim')
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat', 'DataSet')
load('CorrDistribution_Cult6_temp2.mat', 'LagAtMax', 'MaxCosSim')
% MaxCosSim = MaxCosSim{w};
ic= DataSet{w}.ic;
a2a=MaxCosSim(size(ic,2)+1:end,size(ic,2)+1:end);
n2n=MaxCosSim(1:size(ic,2),1:size(ic,2));
a2n=MaxCosSim(1:size(ic,2),size(ic,2)+1:end);
ch = DataSet{w}.channel;
%---Some Top Connections---%
[~,~,~,a1,a_e2,~,ch,~]=ParseConnectionValuesWithVF(w,1000);
astros=[];
for i=1:numel(unique(a_e2))
    [test(i),F(i),~]=mode(a_e2);
    astros = [astros,a1(ismember(a_e2,test(i)))];
    a_e2(a_e2==test(i))=[];
end
[~,~,~,a1,a_e2,~,ch,~]=ParseConnectionValuesWithVF(w,1000);
% [~,~,channels] = FindElecsinVF(ch,MeaMap);
% a_e2 = [a_e2,channels'];
% [~,ac]=nanmax(a2n(:,ismember(ic(1,:),channels)));
% a1 = [a1, ];

test(F<(mean(F)+std(F)))=[];
astros(F<(mean(F)+std(F)))=[];

a2n = a2n(ismember(ic(1,:),a_e2),unique(a1));
n2n = n2n(ismember(ic(1,:),a_e2),ismember(ic(1,:),a_e2));
a2a = a2a(unique(a1),unique(a1));
ic = ic(:,ismember(ic(1,:),a_e2));

MaxCosSim=zeros(size(n2n,1)+size(a2a,1),size(n2n,1)+size(a2a,1));
MaxCosSim(1:size(n2n,1),1:size(n2n,2))=n2n+n2n'+eye(size(n2n));
MaxCosSim(1:size(n2n,1),size(n2n,2)+1:end)=a2n;
MaxCosSim(size(n2n,2)+1:end,1:size(n2n,1))=a2n';
MaxCosSim(size(n2n,1)+1:end,size(n2n,1)+1:end)=a2a+a2a'+eye(size(a2a));

%---ViewField Set of Electrodes---%
[~,~,channels] = FindElecsinVF(ch,MeaMap);
chIdx = find(ismember(ic(1,:),channels));
corrs = arrayfun(@(x) max(MaxCosSim(ismember(ic(1,:),channels),x)), ic(1,:)); %stopped here

%---Random Set of 9---%
[~,~,rand9] = FindElecsinVF(ic(1,randi(size(ic,2))),MeaMap);
rand9Idx = find(ismember(ic(1,:),rand9));

%---Page Rank N2N---%
n2n = n2n+n2n'+eye(size(n2n));
n2n(logical(isnan(n2n)))=0;
% prn = tiedrank(pagerank(n2n));
prn = pagerank(n2n);
prn=prn';

%---Page Rank A2A---%
a2a = a2a+a2a'+eye(size(a2a));
a2a(logical(isnan(a2a)))=0;
% pra = tiedrank(pagerank(a2a));
pra = pagerank(a2a);
pra=pra';

%---Page Rank A2N---%
mat = zeros(size(n2n,1)+size(a2n,2),size(n2n,1)+size(a2n,2));
mat(1:size(ic,2),size(ic,2)+1:end)=MaxCosSim(1:size(ic,2),size(ic,2)+1:end);
mat=mat+mat'+eye(size(mat));
mat(isnan(mat))=0;

% MaxCosSim(isnan(MaxCosSim))=0;
% MaxCosSim = MaxCosSim+MaxCosSim'+eye(size(MaxCosSim));
% pr = tiedrank(pagerank(MaxCosSim));

% pr = tiedrank(pagerank(mat));
pr = pagerank(mat);
prna = pr(1:size(ic,2))';
pran = pr(size(ic,2)+1:end)';

%------%

ingroup = find(ismember(ic(1,:),test));
outgroup = setdiff(1:numel(ic(1,:)),find(ismember(ic(1,:),test)));

ingroupa = find(ismember(unique(a1),unique(astros)));
outgroupa = setdiff(1:numel(unique(a1)),ingroupa);

subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.2 0.05], [0.05 0.05]);
subplot(1,4,1)
title('Correlations');
notBoxPlot([padarray(corrs(ingroup),[0,numel(outgroup)-numel(ingroup)],nan,'post');... %corr of vf elecs to main pairs
            corrs(outgroup)]'); % corr of vf elecs to rest
set(gca,'ylim',[0,1],'XTickLabel', {'VF Elecs to main pairs', 'VF Elecs to rest'});
set(gca,'FontSize',18);
XYrotalabel(45,0);

subplot(1,4,2:4)
title('Centrality');  
notBoxPlot([padarray(prna(chIdx),[0,numel(outgroup)-numel(chIdx)],nan,'post');...       %page rank of vf elecs in a2n
            prna(outgroup);...                                                          %page rank of rest in a2n
%             padarray(prna(rand9Idx),[0,numel(outgroup)-numel(rand9Idx)],nan,'post');... %page rank of random 9 elecs in a2n
            padarray(prna(ingroup),[0,numel(outgroup)-numel(ingroup)],nan,'post');...   %page rank of main elecs in a2n
            padarray(prn(chIdx),[0,numel(outgroup)-numel(chIdx)],nan,'post');...        %page rank of vf elecs in n2n
%             padarray(prn(rand9Idx),[0,numel(outgroup)-numel(rand9Idx)],nan,'post');...  %page rank of random 9 elecs in n2n
            prn(outgroup);...                                                           %page rank of rest of elecs in n2n
            padarray(prn(ingroup),[0,numel(outgroup)-numel(ingroup)],nan,'post');...    %page rank of main elecs in n2n
%             padarray(pra(ingroupa),[0,numel(outgroup)-numel(ingroupa)],nan,'post');...  %page rank of astro's connected to main elecs
%             padarray(pra(outgroupa),[0,numel(outgroup)-numel(outgroupa)],nan,'post');...%page rank of astro's connected to other elecs                                                   
            ]');       
ylabel('PageRank [ranks]');
% set(gca,'XTickLabel', {'A2N VF Elecs','A2N Rest of Elecs', 'A2N Random 9', 'A2N Main Elecs', 'N2N VF Elecs','N2N Random 9','N2N Rest of Elecs','N2N Main Elecs','A2A Main Elecs','A2A Rest of Elecs'});
set(gca,'FontSize',18);
XYrotalabel(45,0);

%---Scatter Ranks---%

% Scores
MaxCosSim = MaxCosSim+MaxCosSim'+eye(size(MaxCosSim));
MaxCosSim(isnan(MaxCosSim))=0;
figure;
subplot(1,2,1)
hold on;
pr = pagerank(MaxCosSim);
pran = pagerank(mat);
scatter(pr(size(ic,2)+1:end)',pagerank(a2a),'b')
scatter(pran(size(ic,2)+1:end),pagerank(a2a),'r')
xlabel('Astrocyte-Neuron');
ylabel('Astrocyte-Astrocyte');
legend({'Complete Corr Matrix','A2N Corr Matrix'})
set(gca,'FontSize',18);
title('Astrocyte');
ylim([0,11].*10^-3);
subplot(1,2,2)
hold on;
scatter(pr(1:size(ic,2))',pagerank(n2n),'k')
scatter(pran(1:size(ic,2))',pagerank(n2n),'g')
xlabel('Astrocyte-Neuron');
ylabel('Neuron-Neuron');
legend({'Complete Corr Matrix','A2N Corr Matrix'})
title('Neuron');
set(gca,'FontSize',18);
xlim([0,7].*10^-3);
ylim([0,11].*10^-3);

% Ranks
pr = tiedrank(pr);
pran = tiedrank(pran);
figure;
subplot(1,2,1)
hold on;
scatter(pr(size(ic,2)+1:end),pra,'b')
scatter(pran(size(ic,2)+1:end),pra,'r')
xlabel('Astrocyte-Neuron Rank');
ylabel('Astrocyte-Astrocyte Rank');
legend({'Complete Corr Matrix','A2N Corr Matrix'})
title('Astrocyte');
set(gca,'FontSize',18);
subplot(1,2,2)
hold on;
scatter(pr(1:size(ic,2))',prn,'k')
scatter(prna',prn,'g')
xlabel('Astrocyte-Neuron Rank');
ylabel('Neuron-Neuron Rank');
legend({'Complete Corr Matrix','A2N Corr Matrix'})
title('Neuron');
set(gca,'FontSize',18);

end