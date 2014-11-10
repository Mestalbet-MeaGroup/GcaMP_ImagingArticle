
%---Figure 3---%
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts.mat');
for 	i=1:size(DataSet,1)
    [numBI(i),numNoB(i),numSB(i),PeaksWithBursts{i},BurstsWithPeaks{i},PeaksNoBursts{i},BurstsNoPeaks{i},PeaksWithSBs{i},SBsWithPeaks{i},peakIndex{i}]=CalcBurstPeakOccurenceRate(DataSet,i);
end
k=6;
pi    = PeaksWithBursts{k}; % Find Set of Peaks with Bursts Nearby
bi    = BurstsWithPeaks{k}; % Find Set of Peaks with Bursts Far Away
Nob_pi = PeaksNoBursts{k}; % Find Set of Peaks with Bursts Far Away
No_bi = BurstsNoPeaks{k}; % Find Set of Bursts with no corresponding Peaks
SB_pi = PeaksWithSBs{k}; % Find Set of Peaks near SuperBursts
SB_bi = SBsWithPeaks{k}; % Find Set of SuperBursts near peaks
%% Subplot 1-2: Astrocyte Trace Raster
PlotAstroTracesWithRaster(DataSet{k}.dfTraces',DataSet{k}.dfTime,DataSet{k}.Trim.t,DataSet{k}.Trim.ic);

%% Subplot 3: Set of Peaks with Bursts Nearby
figure('color','white');
NumSecs = 6;
index = find(cellfun(@(x) numel(x),bi)>3);
whichB=index(randi(numel(index)));

peaks = bi{whichB};
WhichTraces = randperm(length(peaks));
WhichTraces=peaks(WhichTraces(1:3));
start = DataSet{k}.bs(whichB)-2*12000;
stop  = DataSet{k}.bs(whichB)+(NumSecs-2)*12000;

% Find Traces Which Have Peaks around the Same time
% [WhichTraces,start,stop] = CalcBestTrio(pi(1:10),DataSet{k}.dfTraces,DataSet{k}.dfTime,NumSecs,DataSet{k}.fs);
RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,start,stop,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
% export_fig('BurstWithTrace.eps','-r600');
% close all;
%% Subplot 4: Superbursts and Traces
figure('color','white');
NumSecs = 75;
%----Set of Traces to plot----%
whichSB=randi(size(SB_bi,2));
peaks = SB_bi{whichSB};
WhichTraces = randperm(length(peaks));
WhichTraces=peaks(WhichTraces(1:3));
start = DataSet{k}.sbs(whichSB)-5*12000;
stop  = DataSet{k}.sbs(whichSB)+(NumSecs-5)*12000;
RasterPlotLineTrace(DataSet{k}.t,DataSet{k}.ic,start,stop,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
%% Subplot 5: Set of Peaks Far from Bursts
figure('color','white');
NumSecs = 15;
% whichtrace=1;
% start = (DataSet{k}.dfTime(peakIndex{k}{whichtrace}(Nob_pi{whichtrace}(2)))-3)*12000;
% stop = (DataSet{k}.dfTime(peakIndex{k}{whichtrace}(Nob_pi{whichtrace}(2)))+(NumSecs-3))*12000;

Nob_piPos = cellfun(@(x,y) peakIndex{k}{x}(y)',num2cell([1:size(Nob_pi,2)]),Nob_pi,'UniformOutput',false);
[WhichTraces,start,stop] = CalcBestTrio(Nob_piPos(11:30),DataSet{k}.dfTraces,DataSet{k}.dfTime,NumSecs,DataSet{k}.fs);
RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,start,stop,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
%% Subplot 6: Set of Bursts with no Peaks

%some of these still seem to have peaks nearby
figure('color','white');
NumSecs = 6;
WhichBurst = find(cellfun(@(x) numel(x)>2,No_bi));
WhichBurst = WhichBurst(randperm(length(WhichBurst)));
WhichBurst = WhichBurst(1);
WhichTraces = No_bi{WhichBurst};
WhichTraces = WhichTraces(randperm(length(WhichTraces)));
WhichTraces = WhichTraces(1:3);

% start = DataSet{k}.bs(WhichBurst)-1*12000;
% stop= DataSet{k}.bs(WhichBurst)+(NumSecs-1)*12000;

RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,DataSet{k}.bs(WhichBurst)-0.1*NumSecs*12000,DataSet{k}.bs(WhichBurst)+NumSecs*12000,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
%% Subplot 7: Bar chart percentage of occurence
figure('color','white');
hold on;
bars = [[nanmean(numBI);1-nanmean(numBI)],[nanmean(numNoB);1-nanmean(numNoB)],[nanmean(numSB);1-nanmean(numSB)]]';
ebars = repmat([nanstd(numBI),nanstd(numNoB),nanstd(numSB)],2,1);
h= bar(bars,'stacked');
set(h(1),'facecolor',[84,95,255]./255);
set(h(2),'facecolor',[157,118,208]./255);

barsx = get(h,'XData');
barsy = get(h,'YData');

h=errorbar(barsx{1},barsy{1},ebars(1,:),'.r');
% h=errorbar(barsx{2},barsy{2},ebars(:,2),'.');
axis tight;
set(gca,'XTick',[1,2,3],'XTickLabel',{'Bursts With Peaks','Peaks Without Bursts','Superbursts With Peaks'});
set(gca,'FontSize',18,'TickDir','Out');
ylabel('Occurence');
ylim([0,1.01]);

%% Subplot 8: Probability a ROI responds to a network event

