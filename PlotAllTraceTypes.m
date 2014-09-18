% Plot Ca transients in relation to bursts
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat')
for i=1:9
t= DataSet{i}.t;
ic=DataSet{i}.ic;
traces=DataSet{i}.dfTraces;
time=DataSet{i}.dfTime;
bs= DataSet{i}.bs;
be= DataSet{i}.be;
[PeakTypeBurst{i},PeakTypeOther{i},CloseTimes{i},FarTimes{i}]=ClassifyAstroPeaks(t,ic,traces,time,bs);
% [PeakTypeBurst,PeakTypeOther]=ClassifyAstroPeaks(t,ic,traces(1:end/2,:),time(1:end/2),bs,be,i);
% [PeakTypeBurst{i},PeakTypeOther{i}]=GroupTracesNearBurstStart(t,ic,traces,time,bs);
end

%% Plot
for i=1:9
screen_size = get(0, 'ScreenSize');
figure('visible','on','Position', [0 0 screen_size(3) screen_size(4) ])
xlims = [0,60];
ylims = [nanmin([PeakTypeBurst{i}(:);PeakTypeOther{i}(:)]),nanmax([PeakTypeBurst{i}(:);PeakTypeOther{i}(:)])];

temp = zscore(inpaint_nans(PeakTypeBurst{i},1));
if size(temp,2)>size(temp,1)
    g = kmeans(temp',3);
else
    g = ones(size(temp,2),1);
end

subplot(3,2,1); 
plot(PeakTypeBurst{i}(:,g==1)); 
title(['Close ',num2str(nanmean(CloseTimes{i}(g==1)))]);
set(gca,'xLim',xlims,'YLim',ylims)

subplot(3,2,3); 
plot(PeakTypeBurst{i}(:,g==2));
title(num2str(nanmean(CloseTimes{i}(g==2))));
set(gca,'xLim',xlims,'YLim',ylims)

subplot(3,2,5); 
plot(PeakTypeBurst{i}(:,g==3));
title(num2str(nanmean(CloseTimes{i}(g==3))));
set(gca,'xLim',xlims,'YLim',ylims)


temp = zscore(inpaint_nans(PeakTypeOther{i},1));
if size(temp,2)>size(temp,1)
    g = kmeans(temp',3);
else
    g = ones(size(temp,2),1);
end


subplot(3,2,2); 
plot(PeakTypeOther{i}(:,g==1));
title(['Far ',num2str(nanmean(FarTimes{i}(g==1)))]);
set(gca,'xLim',xlims,'YLim',ylims)

subplot(3,2,4); 
plot(PeakTypeOther{i}(:,g==2));
title(num2str(nanmean(FarTimes{i}(g==2))));
set(gca,'xLim',xlims,'YLim',ylims)

subplot(3,2,6); 
plot(PeakTypeOther{i}(:,g==3));
title(num2str(nanmean(FarTimes{i}(g==3))));
set(gca,'xLim',xlims,'YLim',ylims);
% export_fig(['TracesNearBursts_Cult',num2str(i),'.tif'],'-r300');
% close all;
end

