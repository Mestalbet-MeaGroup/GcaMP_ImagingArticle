%---Figure 3---%
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim.mat')
k=6;
%% Subplot 1-2: Astrocyte Trace Raster
PlotAstroTracesWithRaster(DataSet{k}.dfTraces',DataSet{k}.dfTime,DataSet{k}.Trim.t,DataSet{k}.Trim.ic);

%% Find Categories

[peakIndex,~,~,ClosestValues,BurstIndex]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,DataSet{k}.Trim.bs,DataSet{k}.t);

%---Find Set of Bursts and Traces That Co-Occur---%
for i=1:size(peakIndex,2)
dis = ClosestValues{i};
index = find((dis>0)&(dis<1));
bi{i} = BurstIndex{i}(index);
pi{i} = peakIndex{i}(index);
end
[peakIndex1,~,~,ClosestValues1,BurstIndex1]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,sort([DataSet{k}.sbs,DataSet{k}.Trim.bs]),DataSet{k}.t);  %Includes SBs, so peaks within SBs are selected by mistake
%---Find Set of Traces with no Bursts---%
% Fix here. Not finding proper pairs.
sbslocs=find(ismember(sort([DataSet{k}.sbs,DataSet{k}.Trim.bs]),DataSet{k}.sbs)==1); %Superbursts
for i=1:size(peakIndex1,2)
dis = abs(ClosestValues1{i});
% dis(ismember(BurstIndex1{i},sbslocs))=[]; %Don't find peaks within a superburst period
index = find(dis>3);
temp = BurstIndex1{i}(index);
remove = (ismember(temp,sbslocs));
temp(remove)=[];
Nob_bi{i}=temp;
temp = peakIndex1{i}(index);
temp(remove)=[];
Nob_pi{i}=temp;
end


%---Find Set of Traces with no Bursts---%

%---Find Set of Superbursts and Traces That Co-Occur---%
[SB_pi,~,~,SB_cv,SB_bi]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,DataSet{k}.sbs,DataSet{k}.t);
%---Find Set of Bursts with no Traces---%



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
% print('-dpng','-r300','Fig_TraceRasterD.png');
% export_fig('Fig_TraceRaster_D.eps','-r300');
% close all;
%-------------------------%

%% Subplot 4: Superbursts and Traces
%------------SB+Trace-----------%
figure('renderer','zbuffer','visible','off','Position', [0 0 screen_size(3) screen_size(4) ])
NumSecs = 75;
RasterPlotLineTrace(DataSet{k}.t,DataSet{k}.ic,DataSet{k}.sbs(1)-0.1*NumSecs*12000,DataSet{k}.sbs(1)+NumSecs*12000,DataSet{k}.dfTraces(:,15:18)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
% print('-dpng','-r300','Fig_TraceRasterE.png');
export_fig('Fig_TraceRaster_E.eps','-r300');
close all;
%Minute 7
%-------------------------%

%% Subplot 5: No bursts but Traces
% Fix here. Not finding proper pairs.
%------------NoBurst+Trace-----------%
NumSecs = 6;
window = ceil(6*DataSet{k}.fs);
WhichTraces = 27:30;
combs = VChooseK(1:numel(WhichTraces),2);
for i=1:size(combs,1)
    p1 = Nob_pi{WhichTraces(combs(i,1))};
    p2 = Nob_pi{WhichTraces(combs(i,2))};   
    [cv,s2i,s1i]=FindNearestElements(p1',p2');
    p{i} = [p1(s1i(cv<window));p2(s2i(cv<window))];
end
numbers = cellfun(@(x) numel(x),p);    
test=[]; for i=1:6, test=[test;ones(numbers(i),1).*i]; end
test(:,2)=cell2mat(p')';
test = sortrows(test,2);
test(:,3)=abs(diff([test(1,2);test(:,2)]))>window;
[inset,~]=initfin(test(:,3)');
inset= [1,inset];
for i=1:numel(inset)-1
    sets{i} = unique(test(inset(i):inset(i+1),2));
    best(i)=numel(unique(test(inset(i):inset(i+1),1)));
end

start = DataSet{6}.dfTime(min(sets{find(best==max(best),1,'Last')})-ceil(window/2))*12000;
stop = DataSet{6}.dfTime(min(sets{find(best==max(best),1,'Last')})+ceil(window/2))*12000;
RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,start,stop,DataSet{k}.dfTraces(:,WhichTraces)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
%-------------------------%

%% Subplot 6: Bursts but No Traces
%------------Burst+NoTrace-----------%
figure('renderer','zbuffer','visible','off','Position', [0 0 screen_size(3) screen_size(4) ])
NumSecs = 6;
WhichBurst = 10;
RasterPlotLineTrace(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,DataSet{k}.Trim.bs(WhichBurst)-0.1*NumSecs*12000,DataSet{k}.Trim.bs(WhichBurst)+NumSecs*12000,DataSet{k}.dfTraces(:,15:18)',DataSet{k}.dfTime.*12000);
set(gca,'PlotBoxAspectRatio',[1,1,1]);
% print('-dpng','-r300','Fig_TraceRasterC.png');
export_fig('Fig_TraceRaster_C.eps','-r300');
close all;
%Minute 42
%-------------------------%

%% Subplot 7: Bar chart percentage of occurence