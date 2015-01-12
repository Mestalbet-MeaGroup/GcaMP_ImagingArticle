fclose('all');clear all; close all;clc;
% load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts.mat');
% load('BurstPeakOccurance.mat');
% load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts_withSBbursts.mat');
% load('BurstPeakOccurance_withSBbursts.mat');
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts_withSBbursts.mat')
% load('BurstPeakOccurance_withSBbursts_onlyforward.mat');
load('BurstPeakOccurance_withSBbursts_onlyforward_tillNextBurst.mat');

BurstsNoPeaks(4)=[];
BurstsWithPeaks(4)=[];
DataSet(4)=[];
numBI(4)=[];
numNoB(4)=[];
numSB(4)=[];
peakIndex(4)=[];
PeaksNoBursts(4)=[];
PeaksWithBursts(4)=[];
PeaksWithSBs(4)=[];
SBsWithPeaks(4)=[];

load('MeaMapPlot.mat','MeaMap');
MeaMap = rot90(MeaMap,-1);
for i=1:8
    peaks = DataSet{i}.dfTime(cell2mat(peakIndex{i}'));
    validBS{i} = find( (DataSet{i}.bs<=max(peaks)*12000) & (DataSet{i}.bs>=min(peaks)*12000) );
    temp = cellfun(@(x) isempty(x),BurstsWithPeaks{i});
    x{i} = intersect(find(temp),validBS{i});
end

%% Make burst-Ca calculations
c1=1;
c2=1;
noResponseBs=[];
ResponseBs=[];
orderR=[];
orderNR=[];
cultR=[];
cultNR=[];
NBcult=[];
Bcult=[];
for i=1:8;
t= DataSet{i}.t;
ic = DataSet{i}.ic;
bs = DataSet{i}.bs;
be = DataSet{i}.be;
vec=ConvertIC2Samora(ic);

%---Find Bursts with and Without Ca Response---%
y = setdiff(1:numel(bs),x{i});
y = intersect(y,validBS{i});
rB = bs(x{i});
rBe = be(x{i});
nrB = bs(y);
nrBe = be(y);

%---Calculate Distances between Bursts with/without Ca Response---%

% distances = abs(distances);
distances = bsxfun(@minus,rB',[rB,nrB]); %yes response bursts vs. bursts
distances(distances==0)=nan;
% Distances for first preceding burst and first following burst %
dbehind = distances; 
dahead= distances; 
dbehind(dbehind>0)=nan;
dahead(dahead<0)=nan;
[~,noResB]=nanmin(abs(dbehind),[],2); % distance of response bursts to nearest  burst
[~,noResA]=nanmin(abs(dahead),[],2); % distance of response bursts to nearest  burst

for kk=1:numel(noResA); 
    noResA(kk)=dahead(kk,noResA(kk)); 
end

for kk=1:numel(noResB); 
    noResB(kk)=dbehind(kk,noResB(kk)); 
end

% noRes = [noResA;noResB]; %first burst ahead and behind
noRes = [noResA;[]]; %first burst ahead 

% noRes(isnan(noRes))=[];
% distances = bsxfun(@minus,rB',nrB); %no response bursts vs. yes response bursts
% [~,noRes]=min(abs(distances),[],2); % distance of no response bursts to nearest response burst
% for kk=1:numel(noRes); noRes(kk)=distances(kk,noRes(kk)); end

distances = bsxfun(@minus,nrB',[rB,nrB]); %no response bursts vs. bursts
% distances = abs(distances);
% distances = triu(distances,1)+tril(nan(size(distances)));
distances(distances==0)=nan;
dbehind = distances; 
dahead= distances; 
dbehind(dbehind>0)=nan;
dahead(dahead<0)=nan;
[~,ResB]=nanmin(abs(dbehind),[],2); % distance of response bursts to nearest  burst
[~,ResA]=nanmin(abs(dahead),[],2); % distance of response bursts to nearest  burst

for kk=1:numel(ResA); 
    ResA(kk)=dahead(kk,ResA(kk)); 
end

for kk=1:numel(ResB); 
    ResB(kk)=dbehind(kk,ResB(kk)); 
end

% Res = [ResA;ResB]; %first burst ahead and behind
Res = [ResA;[]]; %first burst ahead 

% [~,Res]=nanmin(abs(distances),[],2); % distance of no response bursts to nearest burst
% for kk=1:numel(Res); Res(kk)=distances(kk,Res(kk)); end
% Res(isnan(Res))=[];

noResponseBs = [noResponseBs;noRes]; %summary for all cultures
cultNR = [cultNR;ones(size(noRes)).*i];%for nested anova
ResponseBs = [ResponseBs;Res];
cultR = [cultR;ones(size(Res)).*i];%for nested anova
%---Calculate Spikes+Recruitment (num electrodes) in Responsive Bursts vs. Non-Responsive Bursts---%
[~,~,channels] = FindElecsinVF(DataSet{i}.channel,MeaMap);
numc(i)  = numel(channels);

for j=1:numel(rB);
    frB(c1) = sum(t>=rB(j) & t<=rBe(j))/((rBe(j)-rB(j))/12000)/numel(unique(vec)); %spikes per second in burst
    channelOrderR = unique(vec(t>=(rB(j)*0.9) & t<=(rBe(j)*1.1)),'stable');
    order=[];
    [~,~,order]=intersect(channels,channelOrderR); %spike order which propagates through viewfield for responsive
    orderR = [orderR;order];
    NumVfR{i}(j)=numel(order); %Number of channels in burst which fall within VF
    recruitmentB(c1) = numel(channelOrderR)/numel(unique(vec)); % percentage of channels recruited to burst
    Bcult(c1)=i; %for nested anova
    c1=c1+1;
end

for j=1:numel(nrB);
    frNB(c2) = sum(t>=nrB(j) & t<=nrBe(j))/((nrBe(j)-nrB(j))/12000)/numel(unique(vec)); %spikes per second in burst
    channelOrderNR = unique(vec(t>=(nrB(j)*0.9) & t<=(nrBe(j)*1.1)),'stable');
    order=[];
    [~,~,order]=intersect(channels,channelOrderNR); %spike order which propagates through viewfield for non responsive
    orderNR = [orderNR;order];
    NumVfNR{i}(j)=numel(order); %Number of channels in burst which fall within VF
    recruitmentNB(c2) = numel(channelOrderNR)/numel(unique(vec)); % percentage of channels recruited to burst
    NBcult(c2)=i;%for nested anova
    c2=c2+1;
end

end

frB(recruitmentB==0)=[];
frNB(recruitmentNB==0)=[];

recruitmentNB(recruitmentNB==0)=[];
recruitmentB(recruitmentB==0)=[];

Bcult(recruitmentB==0)=[];
NBcult(recruitmentNB==0)=[];

%% Plots
%---Firing rate of elecs in Burst---%
figure('units','normalized','outerposition',[0 0 1 1])
maxSize = max([numel(frB),numel(frNB)]);
notBoxPlot([padarray(frB',[abs(maxSize-numel(frB)),0],nan,'pre'),padarray(frNB,[0,abs(numel(frNB)-maxSize)],nan,'pre')']);
set(gca,'XTickLabel',{'Bursts w/ Ca response', 'Bursts w/o Ca response'});
ylabel('Spike rate per channel [spikes * sec^-1]');
set(gca,'FontSize',18,'TickDir','Out');
% export_fig('FR_Bursts.eps','-r600');
% close all;

%---Recruitment of elecs in Burst---%
figure('units','normalized','outerposition',[0 0 1 1])
notBoxPlot([padarray(recruitmentB',[abs(maxSize-numel(recruitmentB)),0],nan,'pre'),padarray(recruitmentNB,[0,abs(numel(recruitmentNB)-maxSize)],nan,'pre')']);
set(gca,'XTickLabel',{'Bursts w/ Ca response', 'Bursts w/o Ca response'});
ylabel('Percent active channel participation');
set(gca,'FontSize',18,'TickDir','Out');
% export_fig('Recruitment_Bursts.eps','-r600');
% close all;

%---Time difference between Burst with Ca and nearest burst vs. Burst without Ca and nearest burst ---%
figure('units','normalized','outerposition',[0 0 1 1])
maxSize = max([numel(ResponseBs),numel(noResponseBs)]);
notBoxPlot([padarray(ResponseBs',[0,abs(maxSize-numel(ResponseBs))],nan,'pre')',padarray(noResponseBs,[abs(maxSize-numel(noResponseBs)),0],nan,'pre')]./12000);
set(gca,'XTickLabel',{'Intervals w/ Ca response', 'Intervals w/o Ca response'});
set(gca,'FontSize',18,'TickDir','Out');
ylabel('Latency to nearest burst [s]');
% export_fig('Time_Bursts.eps','-r600');
% close all;

%---Recruitment (percentage) of VF elec's in Burst---%
figure('units','normalized','outerposition',[0 0 1 1])
ratioR=[];ratioNR=[];
for i=1:8
    ratioR = [ratioR,NumVfR{i}/numc(i)];
    ratioNR = [ratioR,NumVfNR{i}/numc(i)];
end

[~,bins]=hist([ratioR,ratioNR],10);
count_numR = histc(ratioR,bins);
count_numR = count_numR./trapz(bins,count_numR);
count_numNR = histc(ratioNR,bins);
count_numNR = count_numNR./trapz(bins,count_numNR);
hold on; 
width1 = 1;
bH1 = bar(bins,count_numR,width1,'b'); 
bH2 = bar(bins,count_numNR,width1/2,'r'); 
xlabel('Percentage view-field participation in bursts');
ylabel('Probability density')
set(gca,'FontSize',18,'TickDir','Out');
% export_fig('VFparticipation_Bursts.eps','-r600');
% close all;

%---Order of VF elec's in Burst---%
figure('units','normalized','outerposition',[0 0 1 1])
[~,bins]=hist([channelOrderR,channelOrderNR],10);
count_coR = histc(channelOrderR,bins);
count_coR = count_coR./trapz(bins,count_coR);
count_coNR = histc(channelOrderNR,bins);
count_coNR = count_coNR./trapz(bins,count_coNR);
hold on; 
width1 = 1;
bH1 = bar(bins,count_coR,width1,'b'); 
bH2 = bar(bins,count_coNR,width1/2,'r'); 
xlabel('Rank orders of firing for view-field electrodes in a burst');
ylabel('Probability density')
legend([bH1,bH2],{'Bursts with Ca','Bursts without Ca'});
set(gca,'FontSize',18,'TickDir','Out');
% export_fig('VForder_Bursts.eps','-r600');
% close all;

% PlotBurstRecruitment(ratioR,ratioNR)
% PlotBurstRecruitmentOrder(channelOrderR,channelOrderNR)
%% Statistics
res = [ResponseBs; noResponseBs]';
gres = [ones(size(ResponseBs));ones(size(noResponseBs)).*2];
cultres=[cultR;cultNR];

fr = [frB,frNB];
rec  = [recruitmentB,recruitmentNB];
g = [ones(size(frB)),ones(size(frNB)).*2];
cult = [Bcult,NBcult];

% clear_all_but('res','gres','cultres', 'fr','rec','g','cult');

for i=1:max(cultres) %latency between bursts
    [p(i),~,~] = ttest2(res(gres==1 & cultres==i),res(gres==2 & cultres==i));
end
[~, ~,pres] = mafdr(p);

for i=1:max(cult)
    [p(i),~,~] = ttest2(rec(g==1 & cult==i),rec(g==2 & cult==i));
end
[~, ~,prec] = mafdr(p);

for i=1:max(cult)
    [p(i),~,~] = ttest2(fr(g==1 & cult==i),fr(g==2 & cult==i));
end
[~, ~,pfr] = mafdr(p);

% nestmat = [ 0,1;0,0];
% [p,table,stats] = anovan(fr,{g,cult},'random',1,'nested',nestmat,'model','full','varnames',{'WithCa' 'Culture'});
% [p,table,stats] = anovan(rec,{g,cult},'random',2,'nested',nestmat,'model','full','varnames',{'WithCa' 'Culture'});
% 
% [p,table,stats] = anovan(res,{gres,cultres},'random',2,'nested',nestmat,'model','full','varnames',{'WithCa' 'Culture'});

%% Probability of a Ca response versus burst characteristics
clear fr; 
for i=1:8
    Pca{i} = cellfun(@(x) numel(x)/size(peakIndex{i},2),BurstsWithPeaks{i});
    if isfield(DataSet{i},'sb_bs')
        bs = sort([DataSet{i}.bs,DataSet{i}.sb_bs]);
        be = sort([DataSet{i}.be,DataSet{i}.sb_be]);
    else
        bs = DataSet{i}.bs;
        be = DataSet{i}.be;
    end
    distances = bsxfun(@minus,bs',bs); 
    distances(distances==0)=nan;
    d = triu(abs(distances),1)+tril(nan(size(distances)));
    latency{i}=min(d,[],2)./12000; 
    
    for j=1:size(BurstsWithPeaks{i},2)
        t=DataSet{i}.t;
        vec=ConvertIC2Samora(DataSet{i}.ic);
        blims = t>=bs(j) & t<=be(j);
        recruitment{i}(j) = numel(unique(vec(blims)));
        fr{i}(j)=sum(blims)/((be(j)-bs(j))/12000)/recruitment{i}(j);
        recruitment{i}(j)= recruitment{i}(j)./numel(unique(vec));
        numbs{i}(j)=sum(abs(distances(j,:))<(5*12000));
    end
    
end
figure;
scatter(cell2mat(fr),cell2mat(Pca));
xlabel('firing rate [Hz]');
ylabel('Probability of an ROI having a Ca peak');
axis tight;

figure;
scatter(cell2mat(recruitment),cell2mat(Pca));
xlabel('recruitment [%]');
ylabel('probability of an ROI having a Ca peak');

figure;
lat = cell2mat(latency');
lat(isnan(lat))=[];
scatter(cell2mat(Pca),lat);
xlabel('latency [s]');
ylabel('probability of an ROI having a Ca peak');

figure;
scatter(cell2mat(numbs),cell2mat(Pca));
xlabel('number bursts within 5 secs [#]');
ylabel('probability of an ROI having a Ca peak');

%% IBI Distributions
plaw=cell(8,1);
distfit=plaw;
dplaw=plaw;
parfor i=1:8
     if isfield(DataSet{i},'sb_bs')
        bs = sort([DataSet{i}.bs,DataSet{i}.sb_bs]);
        be = sort([DataSet{i}.be,DataSet{i}.sb_be]);
    else
        bs = DataSet{i}.bs;
        be = DataSet{i}.be;
     end
    ibi = (be(2:end)-bs(1:end-1))./12000;
%     ibi = diff(ibi);
    dispHIST(ibi,nsOPTBINS(ibi));
if numel(ibi)<300
    plaw{i} = Power_lawFit(ibi,1);
    dplaw{i} = Power_lawFit(diff(ibi),1);
else
    plaw{i} = Power_lawFit(ibi,0);
    dplaw{i} = Power_lawFit(diff(ibi),0);
end
    [distfit{i},~] = allfitdist(ibi);
end
[~, ~,pfr] = mafdr(cellfun(@(x) x.p,plaw));
[~, ~,pfr2] = mafdr(cellfun(@(x) x.p,dplaw));

pooled=[];
for i=1:8
    subplot(4,2,i)
    if isfield(DataSet{i},'sb_bs')
        bs = sort([DataSet{i}.bs,DataSet{i}.sb_bs]);
        be = sort([DataSet{i}.be,DataSet{i}.sb_be]);
    else
        bs = DataSet{i}.bs;
        be = DataSet{i}.be;
    end
    ibi = (be(2:end)-bs(1:end-1))./12000;
    dispHIST(ibi,nsOPTBINS(ibi));
    eval(['ibi',num2str(i),' = (be(2:end)-bs(1:end-1))./12000;']);
    pooled = [pooled,ibi];
end