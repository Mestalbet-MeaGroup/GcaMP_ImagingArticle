
%---Figure 3---%
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS.mat')
k=6;

%------------------------Find Categories of Events----------------------------%
[b2pi,p2bi,p2b_cv,b2p_cv,peaks]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,DataSet{k}.bs,sort(DataSet{k}.t));
[sb2pi,p2sbi,p2sb_cv,sb2p_cv,peaksSB]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,DataSet{k}.sbs,sort(DataSet{k}.t));
% [~,~,p2bSB_cv,~,~]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,[DataSet{k}.bs,DataSet{k}.sbs],sort(DataSet{k}.t));
% b2p = index for each burst of the closest peak
% b2p_cv = corresponding distance between burst and peak

% p2b = index for each peak of the closest burst
% p2b_cv = corresponding distance between peak and burst

% sb2p = index for each superburst of the closest peak
% p2sb_cv = corresponding distance between superburst and peak

% p2sb = index for each peak of the closest superburst
% p2sb_cv = corresponding distance between peak and superburst

%% Find Set of Peaks with Bursts Nearby
for i=1:size(p2bi,2)
    dis = p2b_cv{i};
    index = find((dis>0)&(dis<1));
    bi{i} = p2bi{i}(index);
    pi{i} = peaks{i}(index);
end

%% Find Set of Peaks with Bursts Far Away
for i=1:size(p2bi,2)
    dis = p2b_cv{i};
    index1 = find(abs(dis)>3); 
    dis2 = abs(p2sb_cv{i});
    sbw = DataSet{k}.sbw(p2sbi{i})./12000;
    index2 = find(dis2<sbw);
    index =index1(~ismember(index1,index2));
    Nob_pi{i} = peaks{i}(index);
end

%% Find Set of Bursts with no corresponding Peaks
for i=1:size(p2bi,2)
    dis = b2p_cv{i};
    index1 = find(abs(dis)>3); 
    dis2 = abs(sb2p_cv{i});
    sbw = DataSet{k}.sbw(p2sbi{i}(sb2pi{i}))./12000;
    index2 = find(dis2<sbw);
    index =index1(~ismember(index1,index2));
    No_bi{i} = index;
end
% Are most skipped bursts close to another burst? Are most skipped bursts
% weaker than non-skipped bursts? Farther away?
%% Find Set of Peaks near SuperBursts
for i=1:size(p2sbi,2)
    dis = p2sb_cv{i};
    index = find((dis>0)&(dis<3));
    SB_pi{i} =  peaks{i}(index);
    SB_bi{i} =  unique(p2sbi{i}(index));
end
%% Subplot 1-2: Astrocyte Trace Raster
PlotAstroTracesWithRaster(DataSet{k}.dfTraces',DataSet{k}.dfTime,DataSet{k}.Trim.t,DataSet{k}.Trim.ic);

%% Subplot 3: Set of Peaks with Bursts Nearby
figure;
NumSecs = 6;
%----Set of Traces to plot----%
% WhichTraces = 1:10;
%----Which Peaks to Plot----%
% write code for finding peaks at the same time

%---Select Random Set of Traces around a Random Burst---%
% Find Traces Which Have Peaks around the Same time
[WhichTraces,start,stop] = CalcBestTrio(pi(1:10),DataSet{k}.dfTraces,DataSet{k}.dfTime,NumSecs,DataSet{k}.fs);
RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,start,stop,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
%% Subplot 4: Superbursts and Traces
figure;
NumSecs = 75;
%----Set of Traces to plot----%
% WhichTraces = 1:10;
%----Which Peaks to Plot----%
% write code for finding peaks at the same time

[WhichTraces,start,stop] = CalcBestTrio(SB_pi,DataSet{k}.dfTraces,DataSet{k}.dfTime,NumSecs,DataSet{k}.fs);
RasterPlotLineTrace(DataSet{k}.t,DataSet{k}.ic,start,stop,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
%% Subplot 5: Set of Peaks Far from Bursts
figure;
NumSecs = 15;
%----Set of Traces to plot----%
% WhichTraces = 1:10;
%----Which Peaks to Plot----%
% write code for finding peaks at the same time
[WhichTraces,start,stop] = CalcBestTrio(Nob_pi,DataSet{k}.dfTraces,DataSet{k}.dfTime,NumSecs,DataSet{k}.fs);
RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,start,stop,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
% could be some issues with non detected bursts or strong interburst
% activity.
%% Subplot 6: Set of Bursts with no Peaks
figure;
NumSecs = 6;
[bursts,~,ix]=unique(cell2mat(cellfun(@(x) unique(x),No_bi,'UniformOutput',0)));
n=histc(ix,unique(ix));
[~,WhichBurst]=max(n);
WhichBurst = bursts(ix(WhichBurst));
TracesSet=find(cellfun(@(x) sum(x==WhichBurst)>0,No_bi));
WhichTraces = TracesSet(randperm(numel(TracesSet),3));
RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,DataSet{k}.Trim.bs(WhichBurst)-0.1*NumSecs*12000,DataSet{k}.Trim.bs(WhichBurst)+NumSecs*12000,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
%% Subplot 7: Bar chart percentage of occurence
figure;
numBI = numel(unique(cell2mat(bi)))./numel(DataSet{k}.bs); % Percentage of bursts with a corresponding calcium increase
temp=[];
for i=1:size(Nob_pi,2),
    if ~isempty(Nob_pi{i}),
        temp = [temp;unique(Nob_pi{i})];
    end
end
numNoB = numel(temp)./numel(cell2mat(cellfun(@(x) unique(x),peaks,'UniformOutput',0)')); % Percentage of calcium increases without a corresponding burst
numSB = numel(unique(cell2mat(SB_bi)))./numel(DataSet{k}.sbs); % Percentage of superbursts with a corresponding calcium increase
bar([numBI,numNoB,numSB,1-numBI]);
set(gca,'XTickLabel',{'Bursts-Ca','Ca-No Bursts','SB-Ca','Bursts-No Ca'});
%% Subplot 8: Probability a ROI responds to a network event

