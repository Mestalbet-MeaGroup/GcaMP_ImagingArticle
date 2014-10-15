%---Figure 3---%
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim.mat')
k=6;
%% Subplot 1-2: Astrocyte Trace Raster
PlotAstroTracesWithRaster(DataSet{k}.dfTraces',DataSet{k}.dfTime,DataSet{k}.Trim.t,DataSet{k}.Trim.ic);

%% Find Categories

%---Find Set of Bursts and Traces That Co-Occur---%
[peakIndex,~,~,ClosestValues,BurstIndex]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,DataSet{k}.Trim.bs,DataSet{k}.t);
%--Remove very low amp peaks (avoid type I errors)--%
for i=1:size(DataSet{k}.dfTraces,2)
    amps{i}  = DataSet{k}.dfTraces(peakIndex{i},i);
end
numlevs=3;
amppooled = cell2mat(amps');
ranks = otsu(amppooled,numlevs);
cutoff=max(amppooled(ranks==1));

for i=1:size(peakIndex,2)
    remove=DataSet{k}.dfTraces(peakIndex{i},i)<cutoff;
    peakIndex{i}(remove)=[];
    ClosestValues{i}(remove)=[];
    BurstIndex{i}(remove)=[];
    dis = ClosestValues{i};
    index = find((dis>0)&(dis<1));
    bi{i} = BurstIndex{i}(index); %why are there so few when it looks like it happens all the time?
    pi{i} = peakIndex{i}(index);
end

%---Find Set of Traces with no Bursts---%
[peakIndex1,~,~,ClosestValues1,BurstIndex1]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,sort([DataSet{k}.sbs,DataSet{k}.Trim.bs]),DataSet{k}.t);  %Includes SBs, so peaks within SBs are selected by mistake
%--Remove very low amp peaks (avoid type I errors)--%
for i=1:size(DataSet{k}.dfTraces,2)
    amps{i}  = DataSet{k}.dfTraces(peakIndex1{i},i);
end
numlevs=3;
amppooled = cell2mat(amps');
ranks = otsu(amppooled,numlevs);
cutoff=max(amppooled(ranks==1));

sbslocs=find(ismember(sort([DataSet{k}.sbs,DataSet{k}.Trim.bs]),DataSet{k}.sbs)==1); %Superbursts
for i=1:size(peakIndex1,2)
    remove=DataSet{k}.dfTraces(peakIndex1{i},i)<cutoff;
    peakIndex1{i}(remove)=[];
    ClosestValues1{i}(remove)=[];
    BurstIndex1{i}(remove)=[];
    dis = abs(ClosestValues1{i});
    index = find(dis>3);
    temp = BurstIndex1{i}(index);
    remove = (ismember(temp,sbslocs));
    temp(remove)=[];
    Nob_bi{i}=temp;
    temp = peakIndex1{i}(index);
    temp(remove)=[];
    Nob_pi{i}=temp;
end

%---Find Set of Superbursts and Traces That Co-Occur---%
[SBpeakindex,~,~,SBclosestvalues,SBburstindex]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,DataSet{k}.sbs,DataSet{k}.t);
%--Remove very low amp peaks (avoid type I errors)--%
for i=1:size(DataSet{k}.dfTraces,2)
    amps{i}  = DataSet{k}.dfTraces(SBpeakindex{i},i);
end
numlevs=3;
amppooled = cell2mat(amps');
ranks = otsu(amppooled,numlevs);
cutoff=max(amppooled(ranks==1));
for i=1:size(peakIndex,2)
    remove=DataSet{k}.dfTraces(peakIndex{i},i)<cutoff;
    SBpeakindex{i}(remove)=[];
    SBclosestvalues{i}(remove)=[];
    SBburstindex{i}(remove)=[];
    dis = ClosestValues{i};
    index = find((dis>0)&(dis<1));
    SB_bi{i} = BurstIndex{i}(index);
    SB_pi{i} = SBpeakindex{i}(index);
end

%---Find Set of Bursts with no Traces---%
% for i=1:size(bi,2)
%     temp=1:numel(DataSet{k}.Trim.bs);
%     temp(unique(bi{i}))=[];
%     NoT_bi{i} = temp;
% end
temp=1:numel(DataSet{k}.Trim.bs);
temp(unique(cell2mat(bi)))=[];
NoT_bi = temp;

%% Subplot 3: Bursts and Traces

%------------Burst+Trace-----------%
% figure('renderer','zbuffer','visible','off','Position', [0 0 screen_size(3) screen_size(4) ])
NumSecs = 6;
WhichTraces = 15:18;
[WhichBursts,~,ib]=unique(cell2mat(bi(WhichTraces)));
WhichBursts = WhichBursts(find(arrayfun(@(x) sum(ib==x),unique(ib))==numel(WhichTraces)));
WhichBurst = WhichBursts(randi(numel(WhichBursts),1));
RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,DataSet{k}.Trim.bs(WhichBurst)-0.1*NumSecs*12000,DataSet{k}.Trim.bs(WhichBurst)+NumSecs*12000,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
%-------------------------%

%% Subplot 4: Superbursts and Traces
%------------SB+Trace-----------%
NumSecs = 75;
WhichTraces = 15:18;
[WhichBursts,~,ib]=unique(cell2mat(SB_bi(WhichTraces)));
WhichBurst = WhichBursts(randi(numel(WhichBursts),1));
RasterPlotLineTrace(DataSet{k}.t,DataSet{k}.ic,DataSet{k}.sbs(WhichBurst)-0.1*NumSecs*12000,DataSet{k}.sbs(WhichBurst)+NumSecs*12000,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
%-------------------------%

%% Subplot 5: No bursts but Traces

%------------NoBurst+Trace-----------%
numbers=[];
pkidx=[];
window=86;
traces = zscore(DataSet{k}.dfTraces);
for i=1:110
    if ~isempty(Nob_pi{i})
        numbers=[numbers;ones(numel(Nob_pi{i}),1).*i];
        pkidx = [pkidx;Nob_pi{i}];
    end
end
peaks = [numbers,pkidx];
for i=1:size(peaks,1)
    peaks(i,3)=traces(pkidx(i),numbers(i));
end
peaks = sortrows(peaks,2);

combs = VChooseK(1:size(peaks,1),3);
c1 = combs(:,1);
c2 = combs(:,2);
c3 = combs(:,3);
p1a = peaks(c1,1);
p2a = peaks(c2,1);
p3a = peaks(c3,1);
p1b = peaks(c1,2);
p2b = peaks(c2,2);
p3b = peaks(c3,2);
p1c = peaks(c1,3);
p2c = peaks(c2,3);
p3c = peaks(c3,3);
c1 =zeros(numel(combs(:,1)),1);
c2 =zeros(numel(combs(:,1)),1);
c3 =zeros(numel(combs(:,1)),1);
parfor i=1:numel(c1)
    c1(i) = numel(unique([p1a(i),p2a(i),p3a(i)]));
    temp1  = sort([p1b(i),p2b(i),p3b(i)]);
    c2(i) = temp1(end) - temp1(1);
    c3(i) = median([p1c(i),p2c(i),p3c(i)])/(max([p1c(i),p2c(i),p3c(i)])-min([p1c(i),p2c(i),p3c(i)]));
end

get3 = find(c1==3);
getclose = find(c2<window);
gethigh = find(c3>5);

plot(c2(getclose),c3(getclose),'.')
runIntersect = mintersect(get3,getclose,gethigh);
candidate = find(max(c3(runIntersect))==c3);
WhichTraces= peaks(combs(candidate,:),1);
start = DataSet{k}.dfTime(floor(median(peaks(combs(candidate,:),2)))-ceil(window/2))*12000;
stop = DataSet{k}.dfTime(floor(median(peaks(combs(candidate,:),2)))+ceil(window))*12000;
RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,start,stop,traces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
%-------------------------%


%% Subplot 6: Bursts but No Traces
%------------Burst+NoTrace-----------%
NumSecs = 6;
WhichTraces = 15:18;
WhichBursts = NoT_bi;
WhichBurst = WhichBursts(randi(numel(WhichBursts),1));
RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,DataSet{k}.Trim.bs(WhichBurst)-0.1*NumSecs*12000,DataSet{k}.Trim.bs(WhichBurst)+NumSecs*12000,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);

%-------------------------%

%% Subplot 7: Bar chart percentage of occurence
% pi - bursts and traces which co-occur
% Nob_pi - traces peaks which occur away from bursts
% SB_pi - superbursts and traces which co-occur
% NoT_bi - bursts which occur without traces

bar([numel(cell2mat(pi')),sum(cellfun(@(x) numel(x),Nob_pi)),sum(cellfun(@(x) numel(x),SB_pi)),numel(NoT_bi)]);