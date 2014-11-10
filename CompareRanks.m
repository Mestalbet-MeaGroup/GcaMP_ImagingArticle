load('CorrDistributions2.mat', 'MaxCosSim');
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS.mat');

for i=1:9
    [prn{i},pra{i},prna{i},pran{i}]=CalcPageRanks(MaxCosSim{i},DataSet{i}.ic);
end
% prn = cellfun(@(x) tiedrank(x),prn,'UniformOutput',false);
% prn = cellfun(@(x) x./numel(x),prn,'UniformOutput',false);
prn=cell2mat(prn)';

% pra = cellfun(@(x) tiedrank(x),pra,'UniformOutput',false);
% pra = cellfun(@(x) x./numel(x),pra,'UniformOutput',false);
pra=cell2mat(pra)';

% prna = cellfun(@(x) tiedrank(x),prna,'UniformOutput',false);
% prna = cellfun(@(x) x./numel(x),prna,'UniformOutput',false);
prna=cell2mat(prna)';

% pran = cellfun(@(x) tiedrank(x),pran,'UniformOutput',false);
% pran = cellfun(@(x) x./numel(x),pran,'UniformOutput',false);
pran=cell2mat(pran)';
subplot = @(m,n,p) subtightplot (m, n, p, [0.08 0.08], [0.05 0.05], [0.08 0.08]);
subplot(2,2,1);
[counts,bins]=CalcHist(pra);
bar(bins,counts,'histc'); 
title('Ranks in A2A');
set(gca,'FontSize',18);
subplot(2,2,2);
[counts,bins]=CalcHist(prn);
bar(bins,counts,'histc'); 
title('Ranks in N2N');
set(gca,'FontSize',18);
subplot(2,2,3)
[counts,bins]=CalcHist(pran);
bar(bins,counts,'histc'); 
title('Astro Ranks in A2N');
set(gca,'FontSize',18);
subplot(2,2,4)
[counts,bins]=CalcHist(prna);
bar(bins,counts,'histc'); 
title('Neuron Ranks in A2N');
set(gca,'FontSize',18);
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.005], [0.05 0.05], [0.05 0.05]);

[~,bins] = hist3([prn,prna],[40,50]);
subplot(1,2,1);
hist3([prn,prna],'edges',bins);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
axis tight;
view([0,90]);
ylim([min(bins{2}),max(bins{2})]);
xlim([min(bins{1}),max(bins{1})]);
set(gca,'FontSize',18);

subplot(1,2,2);
hist3([pra,pran],'edges',bins);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
axis tight;
view([0,90]);
set(gca,'YTick',[]);
ylim([min(bins{2}),max(bins{2})]);
xlim([min(bins{1}),max(bins{1})]);
set(gca,'FontSize',18);


steps=0.05;
window=0.1;
nns=[]; Enns=[]; asns=[]; Easns=[]; aas=[]; Eaas=[]; nas=[]; Enas=[];
for i=0:floor(max(prn)/steps)

    start  = 0+i*steps;
    stop = start+window;
    nns(i+1) = mean( prn( (prn>=start) & (prn<=stop) ) );
    Enns(i+1) = std( prn( (prn>=start) & (prn<=stop) ) );
    asns(i+1) = mean( prna( (prn>=start) & (prn<=stop) ) );
    Easns(i+1) = std( prna( (prn>=start) & (prn<=stop) ) );
end

for i=0:floor(max(pra)/steps)
    start  = 0+i*steps;
    stop = start+window;
    aas(i+1) = mean( pra( (pra>=start) & (pra<=stop) ) );
    Eaas(i+1) = std( pra( (pra>=start) & (pra<=stop) ) );
    nas(i+1) = mean( pran( (pra>=start) & (pra<=stop) ) );
    Enas(i+1) = mean( pran( (pra>=start) & (pra<=stop) ) );
end
hold on;
h=plotXYerrorbars(nns, asns, Enns, Easns,'r');
h=plotXYerrorbars(aas, nas, Eaas, Enas,'k');
