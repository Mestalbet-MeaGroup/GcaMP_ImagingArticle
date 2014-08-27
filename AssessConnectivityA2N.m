fclose('all');clear all; close all;clc;
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
set(0,'RecursionLimit',1000);
f = savgol(30,5,0);

%--Mutual Information---%
for i=1:size(DataSet,1)
    traces=DataSet{i}.dfTraces;
    if size(traces,1)<size(traces,2)
        traces=traces';
    end
    tr=filtfilt(f,1,traces);
    fr= DataSet{i}.FR;
    parfor j=1:size(tr,2)
        score(j)= CalcMI(tr(:,j),fr);
    end
    clust = kmeans(score',2);
    results(i)=max([mean(score(clust==1)),mean(score(clust==2))]);
end

%--Correlation (with Entropy filter)---%
for i=1:size(DataSet,1)
    J = entropyfilt(DataSet{i}.dfTraces);
    [tr,b]=CalcDf_f(J',DataSet{i}.fs,DataSet{i}.dfTime);
    tr=tr';
    gfr=histc(sort(DataSet{i}.t)./12000,b);
    parfor j=1:size(tr,2)
        [~,b]=CalcCorrRaw(zscore(tr(:,j)),zscore(gfr),5000);
        score(j)=max(b);
    end
    clust = kmeans(score',2);
    results(i)=max([mean(score(clust==1)),mean(score(clust==2))]);
end

%----Correlation----%
for i=1:size(DataSet,1)
    traces=DataSet{i}.dfTraces;
    if size(traces,1)<size(traces,2)
        traces=traces';
    end
    tr=filtfilt(f,1,traces);
    parfor j=1:size(tr,2)
        [~,b]=CalcCorrRaw(zscore(tr(:,j)),zscore(DataSet{i}.GFR),5000);
        score(j)=max(b);
    end
    clust = kmeans(score',2);
    results(i)=max([mean(score(clust==1)),mean(score(clust==2))]);
end
% CorrGroups = kmeans(results',2);
% if mean(results(CorrGroups==1))>mean(results(CorrGroups==2))
%     CorrGroups(CorrGroups==1)=0;
%     CorrGroups(CorrGroups==2)=1;
%     CorrGroups(CorrGroups==0)=2;
% end
CorrGroups=results;
clear_all_but('CorrGroups','DataSet');
%----Coherence----%
for i=1:size(DataSet,1)
    traces=DataSet{i}.dfTraces;
    time = DataSet{i}.dfTime;
    gfr = DataSet{i}.GFR;
    if size(traces,1)<size(traces,2)
        traces=traces';
    end
    f = savgol(30,3,0);
    tr=filtfilt(f,1,traces);
    score =max(mscohere(zscore(gfr),zscore(tr),[],[],[],DataSet{i}.fs),[],1);
    clust = kmeans(score',2);
    results(i)=max([mean(score(clust==1)),mean(score(clust==2))]);
end
% CohGroups = kmeans(results',2);
% if mean(results(CohGroups==1))>mean(results(CohGroups==2))
%     CohGroups(CohGroups==1)=0;
%     CohGroups(CohGroups==2)=1;
%     CohGroups(CohGroups==0)=2;
% end
CohGroups=results;
clear_all_but('CorrGroups','CohGroups','DataSet');
%----LinearPred----%
for i=1:size(DataSet,1)
    traces=DataSet{i}.dfTraces;
    gfr = DataSet{i}.GFR;
    if size(traces,1)<size(traces,2)
        traces=traces';
    end
    f = savgol(30,3,0);
    tr=filtfilt(f,1,traces);
    parfor j=1:size(tr,2)
        score(j) = LinearPredictionA2N(zscore(gfr),zscore(tr(:,j)));
    end
    clust = kmeans(score',2);
    results(i)=max([mean(score(clust==1)),mean(score(clust==2))]);
end
% LinPredGroups = kmeans(results',2);
% if mean(results(LinPredGroups==1))>mean(results(LinPredGroups==2))
%     LinPredGroups(LinPredGroups==1)=0;
%     LinPredGroups(LinPredGroups==2)=1;
%     LinPredGroups(LinPredGroups==0)=2;
% end
LinPredGroups = results;
clear_all_but('CorrGroups','CohGroups','LinPredGroups','DataSet');
%---Burst Proximity---%
for i=1:size(DataSet,1)
    traces=DataSet{i}.dfTraces;
    time=DataSet{i}.dfTime;
    bs=DataSet{i}.bs;
    [~,~,~,ClosestValues]=FindPeaksNearestBursts(traces,time,bs./12000);
    score(i)=nanmean(ClosestValues);
end
% BurstGroups = kmeans(score',2);
% if mean(score(BurstGroups==1))>mean(score(BurstGroups==2))
%     BurstGroups(BurstGroups==1)=0;
%     BurstGroups(BurstGroups==2)=1;
%     BurstGroups(BurstGroups==0)=2;
% end
BurstGroups = score;
clear_all_but('CorrGroups','CohGroups','LinPredGroups','BurstGroups ','DataSet');
RealGroups3 = [1     1     2     1     2     3     3     2     3]';
RealGroups2 = [1     1     2     1     2     2     2     2     2]';

b = bar([CorrGroups./max(CorrGroups);CohGroups./max(CohGroups);LinPredGroups./max(LinPredGroups);BurstGroups./max(BurstGroups)]'); 
legend(b,{'Correlation','Coherence','Linear Prediction','Burst Proximity'});
% 
% DIV = cellfun(@(x) (x.DIV), DataSet);
% bar(DIV);
% 
% [~,CorrRank]=sort(CorrGroups,'descend');
% [~,CohRank]=sort(CohGroups,'descend');
% [~,LinRank]=sort(LinPredGroups,'descend');
% [~,BurstRank]=sort(BurstGroups,'descend');
% b = bar([CorrRank;CohRank;LinRank;BurstRank]); 
% legend(b,{'Correlation','Coherence','Linear Prediction','Burst Proximity'});