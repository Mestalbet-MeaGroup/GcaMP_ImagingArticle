%% Load Data
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS.mat');
k=6;

%%
[~,peakIndex,~]=CalcPeakStartEnd(DataSet{k}.dfTraces);
% peaks = DataSet{k}.dfTime(unique(cell2mat(peakIndex')));
% whichTrace = cellfun(@(x) numel(x),peakIndex);


%%
counter=1;cf=1;
for i=1:size(peakIndex,2)
    peaks = DataSet{k}.dfTime(peakIndex{i});
    distances = bsxfun(@minus,peaks',DataSet{k}.bs./12000);
    ampvDist{i} = min(abs(distances),[],2);
    tr  = zscore(DataSet{k}.dfTraces(:,i));
    ampvDist{i}(:,2)=tr(peakIndex{i});
    
    
    
    closeby = zeros(size(distances));
    closeby = abs(distances)<1;
    far = zeros(size(distances));
    far = abs(distances)<=3;
    [pc{i},~]=find(closeby==1);
    pf{i}=find(sum(far==0,2)==size(distances,2));
    
    if ~isempty(pc{i})
        pc{i}=unique(pc{i});
        for j=1:numel(pc{i})
            CutoutNear(counter,:) = CreateTraceCutOuts(DataSet{k}.dfTraces(:,i),peakIndex{i}(pc{i}(j)),50,50);
            counter=counter+1;
        end
    end
    
    if ~isempty(pf{i})
        pf{i}=unique(pf{i});
        for j=1:numel(pf{i})
            CutoutFar(cf,:) = CreateTraceCutOuts(DataSet{k}.dfTraces(:,i),peakIndex{i}(pf{i}(j)),50,50);
            cf=cf+1;
        end
    end
    
    
end
ampvDist = cell2mat(ampvDist');
ampvDist =sortrows(ampvDist,1);
%% Subplot 1: Peaks near bursts by amplitude

figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.02], [0.05 0.05], [0.05 0.05]);

numlevs=3;
ampsC  = nanmax(CutoutNear,[],2);
ampsF  = nanmax(CutoutFar,[],2);
ranksC = otsu(ampsC,numlevs);
ranksF = otsu(ampsF,numlevs);
for i=1:numlevs-1
    subplot(numlevs-1,2,2*i)
    plot(CutoutNear(ranksC==i+1,:)','.-');
    ylim([0,max(ampsC)]);
    xlim([0,100]);
    if i==1
        title('Peaks near Bursts')
        set(gca,'XTick',[]);  
    end
    set(gca,'YTick',[]);
    set(gca,'FontSize',18);
    
    subplot(numlevs-1,2,2*i-1)
    plot(CutoutFar(ranksF==i+1,:)','.-');
    ylim([0,max(ampsC)]);
    xlim([0,100]);
    if i==1
        title('Peaks Far from Bursts')
        set(gca,'XTick',[]);
    end
    if i>2
        set(gca,'YTick',[]);
    end
    set(gca,'FontSize',18);
end

%% Subplot 2: Amplitude versus offset
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.05 0.05], [0.05 0.05]);
subplot(6,6,2:6)
hist(ampvDist(:,1),100);
xlim([0,35]);
view([0,90]);
axis off;
set(gca,'YTick',[],'XTick',[]);
set(allchild(gca),'FaceColor',[0 0 0]);
set(gca,'FontSize',18);

subplot(6,6,[7,13,19,25,31])
hist(ampvDist(:,2),100);
axis tight;
xlim([0.5,18]);
set(gca,'YTick',[],'YDir','reverse','XDir','reverse');
box off;
view([90,90]);
set(allchild(gca),'FaceColor',[0 0 0]);
set(gca,'FontSize',18,'ycolor','w');

subplot(6,6,[8:12,14:18,20:24,26:30,32:36])
hist3(ampvDist,[100,100]);
az=0;%35
el=90;%32
view([az,el]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
xlim([0,35]);
ylim([0.5,18]);
set(gca,'YTick',[]);
set(gca,'FontSize',18);

%-----with log scale-----%
hist3(ampvDist,'edges',{logspace(0+eps,max(ampvDist(:,2)),100),linspace(0,max(ampvDist(:,2)),100)});
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
az=0;%35
el=90;%32
view([az,el]);
set(gca,'YTick',[],'XScale','log');
xlim([1,100]);