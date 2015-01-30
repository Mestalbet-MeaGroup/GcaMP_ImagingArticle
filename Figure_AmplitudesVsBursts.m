%% Load Data
% fclose('all');clear all; close all;clc;
% load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts.mat')
% load('BurstPeakOccurance.mat')
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts_withSBbursts.mat')
% load('BurstPeakOccurance_withSBbursts_onlyforward.mat');
load('BurstPeakOccurance_withSBbursts_onlyforward_tillNextBurst.mat');
%---Remove Culture 4 which has virtually no transients---%
BurstsNoPeaks(4)=[];
BurstsWithPeaks(4)=[];
DataSet(4)=[];
numBI(4)=[];
numNoB(4)=[];
numSB(4)=[];
peakIndex(4)=[]; %for each
PeaksNoBursts(4)=[];
PeaksWithBursts(4)=[];
PeaksWithSBs(4)=[];
SBsWithPeaks(4)=[];
%%
AvD = [];AvDcultNum=[];
counter=1;cf=1;
AreaAmp1=[];AreaAmp2=[];CutoutNear=[];CutoutFar=[];NearROIs=[];FarROIs=[];
for k=1:size(DataSet,1)
    %     [~,peakIndex{k},~]=CalcPeakStartEnd(DataSet{k}.dfTraces); %calc instead of load
    %     [~,~,~,PeaksWithBursts{k},~,PeaksNoBursts{k},~,~,~,peakIndex{k}]=CalcBurstPeakOccurenceRate(DataSet,k);
    pi  = peakIndex{k};
    pb  = PeaksWithBursts{k};
    pnb = PeaksNoBursts{k};
    if isfield(DataSet{k},'sb_bs')
        DataSet{k}.bs=sort([DataSet{k}.bs,DataSet{k}.sb_bs]);
    end
    
    for i=1:size(pi,2)
        peaks = DataSet{k}.dfTime(pi{i});
        distances = bsxfun(@minus,peaks',DataSet{k}.bs./12000);
        distances(distances<0)=nan; %look only for peaks which occur after bursts
        ampvDist{i} = nanmin(abs(distances),[],2);
        tr  = DataSet{k}.dfTraces(:,i)-min(DataSet{k}.dfTraces(:,i));
        trz  = zscore(tr);
        ampvDist{i}(:,2)=trz(pi{i});
        Close  = find(ampvDist{i}(:,1)<3);
        far    = find(ampvDist{i}(:,1)>10);
        if ~isempty(Close)
            for j=1:numel(Close)
                CutoutNear = [CutoutNear;CreateTraceCutOuts(tr,pi{i}(Close(j)),100,100)'];
                NearROIs   =  [NearROIs;[k,i,Close(j)]];
            end
        end
        
        if ~isempty(far)
            for j=1:numel(far)
                CutoutFar = [CutoutFar;CreateTraceCutOuts(tr,pi{i}(far(j)),100,100)'];
                FarROIs   =  [FarROIs;[k,i,far(j)]];
            end
        end
        %         if ~isempty(pb{i})
        %             pb{i}=unique(pb{i});
        %             for j=1:numel(pb{i})
        %
        %                 CutoutNear(counter,:) = CreateTraceCutOuts(trz,pi{i}(pb{i}(j)),50,50);
        %                 coNear(counter,:)=[k,i,pb{i}(j)];
        %                 NearROIs(counter,:) = [k,i,pb{i}(j)];
        %                 counter=counter+1;
        %             end
        %         end
        %
        %         if ~isempty(pnb{i})
        %             pnb{i}=unique(pnb{i});
        %             for j=1:numel(pnb{i})
        %                 CutoutFar(cf,:) = CreateTraceCutOuts(trz,pi{i}(pnb{i}(j)),50,50);
        %                 FarROIs(cf,:) = [k,i,pnb{i}(j)];
        %                 cf=cf+1;
        %             end
        %         end
        % for area vs. num pixels
        AreaAmp2 = [AreaAmp2, max(trz)];
    end
    c = regionprops(~DataSet{k}.mask,'Area');
    AreaAmp1 = [AreaAmp1,[c.Area]];
    
    
    temp=cell2mat(ampvDist');
    AvD = [AvD;sortrows(temp,1)];
    AvDcultNum = [AvDcultNum;ones(size(temp,1),1).*k];
end
%% Find ROIs with only Far Peaks or Only Near Peaks or Both
OnlyNearAreas=[];
OnlyFarAreas=[];
NearFarAreas=[];
for r=1:size(DataSet,1)
    [NearFar{r},OnlyNear{r},OnlyFar{r}]=CalcAmpRegions(NearROIs,FarROIs,r);
    c = regionprops(~DataSet{r}.mask,'PixelIdxList','Area');
    area = [c.Area];
    image=double(~DataSet{r}.mask);
    
    NearCounts(r)    = numel(OnlyNear{r});
    FarCounts(r)     = numel(OnlyFar{r});
    NearFarCounts(r) = numel(NearFar{r});
    
    for i=1:numel(OnlyNear{r})
        image(c(OnlyNear{r}(i)).PixelIdxList)=2;
    end
    
    for i=1:numel(OnlyFar{r})
        image(c(OnlyFar{r}(i)).PixelIdxList)=3;
    end
    
    for i=1:numel(NearFar{r})
        image(c(NearFar{r}(i)).PixelIdxList)=4;
    end
    colors =[0,0,0; 1,1,1; 1,0,0; 0,0,1; 0,1,0];
    
    figure;
    %     s1 = subplot(1,2,1);
    imshow(image+1,colors);
    set(gca,'YDir','normal');
    freezeColors;
    
    set(gca,'YDir','normal');
    %     figure;
    OnlyNearAreas = [OnlyNearAreas,area(OnlyNear{r})];
    OnlyFarAreas =  [OnlyFarAreas,area(OnlyFar{r})];
    NearFarAreas =  [NearFarAreas,area(NearFar{r})];
end

%% Plot Areas and Counts versus Ca oscillation
figure;
subplot(2,1,1);
maxsize=max([numel(OnlyNearAreas),numel(OnlyFarAreas),numel(NearFarAreas)]);
% notBoxPlot([padarray(OnlyNearAreas,[0,maxsize-numel(OnlyNearAreas)],nan,'pre');padarray(OnlyFarAreas,[0,maxsize-numel(OnlyFarAreas)],nan,'pre');padarray(NearFarAreas,[0,maxsize-numel(NearFarAreas)],nan,'pre')]');
areas = [padarray(OnlyNearAreas,[0,maxsize-numel(OnlyNearAreas)],nan,'pre');padarray(OnlyFarAreas,[0,maxsize-numel(OnlyFarAreas)],nan,'pre');padarray(NearFarAreas,[0,maxsize-numel(NearFarAreas)],nan,'pre')]';
MeanAreas=nanmean(areas,1);
StdAreas = nanstd(areas,[],1);
hold on;
bar(MeanAreas);
errorbar([1,2,3],MeanAreas,StdAreas,'.k','LineWidth',2);
axis tight;
set(gca,'XTick',[],'TickDir','out');
ylabel('number of pixels per ROI');
subplot(2,1,2);
denominator = NearCounts+FarCounts+NearFarCounts;
counts = [NearCounts;FarCounts;NearFarCounts]./repmat(denominator,3,1);
MeanCounts = nanmean(counts,2);
StdCounts= nanstd(counts,[],2);
hold on;
bar(MeanCounts);
errorbar([1,2,3],MeanCounts,StdCounts,'.k','LineWidth',2);
axis tight;
set(gca,'XTick',[1,2,3],'XTickLabel',{'Ca-burst contemporaneous','Ca-burst non-contemporaneous','Ca-burst both'},'TickDir','out');
ylabel('ratio of ROIs')

%---Stat Tests---%
% a1 = areas(~isnan(areas(:,1)),1);
% a2 = areas(~isnan(areas(:,2)),1);
% a3 = areas(~isnan(areas(:,3)),1);
% a=[a1;a2;a3];
% g1=[ones(size(a1));ones(size(a2))*2;ones(size(a3)).*3];
% [~,~,stats]=anovan(a,g1);
% table = multcompare(stats);
% c1 = counts(1,:)';
% c2 = counts(2,:)';
% c3 = counts(3,:)';
% c=[c1;c2;c3];
% g2=[ones(size(c1));ones(size(c2))*2;ones(size(c3)).*3];
% [~,~,stats]=anovan(c,g2);
% table = multcompare(stats)
%% 2D histogram of amplitude versus region size

figure;
% optM = nsOPTBINS([AreaAmp1;AreaAmp2]);
[histo,bins]=hist3([AreaAmp1;AreaAmp2]','edges',{linspace(0,max(AreaAmp1),20),linspace(0,max(AreaAmp2),20)});
normhist = histo./trapz(bins{2},trapz(bins{1},histo));
xo = repmat(bins{1}',1,size(normhist,2));
yo = repmat(bins{2},size(normhist,1),1);
s = surface(xo,yo,normhist);
axis tight
az=0;%35
el=90;%32
view([az,el]);
set(gca,'PlotBoxAspectRatio',[numel(bins{1}),numel(bins{2}),1],'TickDir','out');
xlabel('number of pixels','FontSize',9);
ylabel('relative Ca amplitude [z-score]','FontSize',9);
c = colorbar;
ylabel(c,'probability density','FontSize',9);
set(c,'TickDir','out');


%% Peaks near bursts by amplitude

figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.02], [0.05 0.05], [0.05 0.05]);
numlevs=4;
ampsC  = nanmax(CutoutNear,[],2);
ampsF  = nanmax(CutoutFar,[],2);

subplot(1,2,1)
hold on
ClusterIDs = otsu(ampsC,numlevs);
color = cool(max(ClusterIDs)+5);
for i=1:max(ClusterIDs)
    plotmat = mean(CutoutNear(ClusterIDs==i,:)',2);
    plot(plotmat,'.-','Color',color(i,:));   
end
axis tight;
ylimits = get(gca,'YLim');

subplot(1,2,2)
hold on
ClusterIDs = otsu(ampsF,numlevs);
color = cool(max(ClusterIDs)+5);
for i=1:max(ClusterIDs)
    plotmat = mean(CutoutFar(ClusterIDs==i,:)',2);
    plot(plotmat,'.-','Color',color(i,:));   
end
axis tight;
ylim(ylimits);

ranksC = otsu(ampsC,numlevs);
ranksF = otsu(ampsF,numlevs);
for i=1:numlevs-1
    subplot(numlevs-1,2,2*i)
    hold on;
    plot((CutoutNear(ranksC==i+1,:))','.-');
    p=plot(mean(CutoutNear(ranksC==i+1,:))','--k','LineWidth',4);
    uistack(p,'top');
    ylim([0,max(ampsC)]);
%     xlim([0,100]);
    if i==1
        title('Peaks near Bursts')
        set(gca,'XTick',[]);
    end
    set(gca,'YTick',[]);
    set(gca,'FontSize',18);
    
    subplot(numlevs-1,2,2*i-1)
    hold on
    plot((CutoutFar(ranksF==i+1,:))','.-');
    p = plot(mean(CutoutFar(ranksF==i+1,:))','--k','LineWidth',4);
    uistack(p,'top');
    ylim([0,max(ampsC)]);
%     xlim([0,100]);
    if i==1
        title('Peaks Far from Bursts')
        set(gca,'XTick',[]);
    end
    if i>2
        set(gca,'YTick',[]);
    end
    set(gca,'FontSize',18);
end

%% Amplitude versus offset
close all; clear bins;
AvDcultNum(AvD(:,1)>25,:)=[];
AvD(AvD(:,1)>25,:)=[];
% AvD(AvD(:,2)<1,:)=[];
figure('Color','white');
% subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.05 0.05], [0.05 0.05]);
[nm,~] = PlotColWiseNormHist(AvD(:,1),AvD(:,2),[100,100],'both');
% for i=1:8
%     [normhist(:,:,i),bins] = PlotColWiseNormHist(AvD(AvDcultNum==i,1),AvD(AvDcultNum==i,2),[100,100],'both');
%     close all;
% end
% xo = repmat(bins{1}',1,size(normhist,2));
% yo = repmat(bins{2},size(normhist,1),1);
% surface(xo,yo,nanmean(normhist,3));
% az=0;%35
% el=90;%32
% view([az,el]);
% set(get(gca,'child'),'FaceColor','flat','CDataMode','auto');
% set(gca,'YTick',[],'TickDir','out','FontSize',9);
% colormap([[0,0,0];parula(numel(unique(normhist(:))))]);
set(findall(gcf,'-property','FontSize'),'FontSize',18);


%--Difference over expected value---%
figure('Color','white');
optM = [100,100];
PlotHistOverExpected(AvD,optM);
xlabel('time from burst start [s]');
ylabel('probability over expected');
set(findall(gcf,'-property','FontSize'),'FontSize',18);

%% Amp vs Num bursts
c=1;
win = 10;
for i=1:size(DataSet,1)
    if isfield(DataSet{i},'sb_bs')
        bs = sort([DataSet{i}.bs,DataSet{i}.sb_bs])./12000;
        be = sort([DataSet{i}.be,DataSet{i}.sb_be])./12000;
    else
        bs = DataSet{i}.bs./12000;
        be = DataSet{i}.be./12000;
    end
    peakTimes{i} = cellfun(@(x) DataSet{i}.dfTime(x),peakIndex{i},'UniformOutput',false);
    for j=1:size(peakTimes{i},2)
        trz  = zscore(DataSet{i}.dfTraces(:,j)-min(DataSet{i}.dfTraces(:,j)));
        for k=1:numel(peakIndex{i}{j})
            AmpVnumB(c,1) = trz(peakIndex{i}{j}(k)); %peak amplitude
            AmpVnumB(c,2) = sum((bs>=(peakTimes{i}{j}(k)-win)) & (bs<=(peakTimes{i}{j}(k)+win))); %number of bursts within window
            c=c+1;
        end
    end
end
close all;
figure;
PlotColWiseNormHist(AmpVnumB(:,2),AmpVnumB(:,1),[20,20],'both');
% [temphist,bins]=hist3(AmpVnumB(:,[2,1]),[20,100]);
% normhist = temphist/trapz(bins{2},trapz(bins{1},temphist));
% xo = repmat(bins{1}',1,size(normhist,2));
% yo = repmat(bins{2},size(normhist,1),1);
% surface(xo,yo,normhist);
% az=0;%35
% el=90;%32
% view([az,el]);
% set(get(gca,'child'),'FaceColor','flat','CDataMode','auto');
% axis tight;
ylabel('peak amplitudes');
xlabel(['# bursts within ',num2str(win), ' s of peak']);