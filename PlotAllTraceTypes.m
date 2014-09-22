% Plot Ca transients in relation to bursts
clear all;
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat')
PeakTypeBurst=cell(9,1);
PeakTypeOther=PeakTypeBurst;
CloseTimes=PeakTypeBurst;
FarTimes=PeakTypeBurst;
parfor i=1:9
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
    xlims = [0,100];
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
%% Pooled PDF
test=[]; for i=1:9; test=cat(2,test,PeakTypeBurst{i}); end
temp = zscore(inpaint_nans(test,1));
ampsNear = nanmax(test,[],1);
g = kmeans(test',3);

subplot(3,2,1);
plot(test(:,g==1));

subplot(3,2,3);
plot(test(:,g==2));

subplot(3,2,5);
plot(test(:,g==3));

test=[]; for i=1:9; test=cat(2,test,PeakTypeOther{i}); end
temp = zscore(inpaint_nans(test,1));
ampsFar = nanmax(test,[],1);
g = kmeans(test',3);

subplot(3,2,2);
plot(test(:,g==1));

subplot(3,2,4);
plot(test(:,g==2));

subplot(3,2,6);
plot(test(:,g==3));

figure;
s1=subplot(2,1,1);
s2=subplot(2,1,2);

h = dispHIST(ampsNear, OPTBINS(ampsNear,1000));
copyobj(allchild(gca),s1);
close(h);
h = dispHIST(ampsFar, OPTBINS(ampsFar,1000));
copyobj(allchild(gca),s2);
close(h);
title(s1,'Amplitudes Near Bursts');
title(s2,'Amplitudes Far from Bursts');
ylabel(s1,'Mean Bin Probabilities');
ylabel(s2,'Mean Bin Probabilities');
xlabel(s2,'Amplitudes');
set(s1,'XLim',[0,4],'YLim',[-0.1,8]);
set(s2,'XLim',[0,4],'YLim',[-0.1,8]);
%% Pooled CDF

hold on;
[CdfY,CdfX,flo,fup] = ecdf(ampsNear','Function','cdf');
CdfX=CdfX(2:end)';
CdfY=CdfY(2:end)';
flo=flo(2:end)';
fup=fup(2:end)';
X = [CdfX,fliplr(CdfX)];
Y = [flo(1:end-1),CdfY(end),CdfY(end), fliplr(fup(1:end-1))];
plot(CdfX,CdfY,'ok','MarkerSize',4,'MarkerFaceColor','k');
patch(X,Y,'b','FaceAlpha',0.2,'EdgeColor','b');
xPred = fitdist(SBibi', 'binomial');
[yPred,YLower,YUpper] = cdf(xPred,CdfX);
plot(CdfX,yPred,'-r','LineWidth',3);

[CdfY,CdfX,flo,fup] = ecdf(ampsFar','Function','cdf');
CdfX=CdfX(2:end)';
CdfY=CdfY(2:end)';
flo=flo(2:end)';
fup=fup(2:end)';
X = [CdfX,fliplr(CdfX)];
Y = [flo(1:end-1),CdfY(end),CdfY(end), fliplr(fup(1:end-1))];
plot(CdfX,CdfY,'sk','MarkerSize',4,'MarkerFaceColor','k');
patch(X,Y,'r','FaceAlpha',0.2,'EdgeColor','r');

legend({'','Near','','Far'});

%% Seperate CDF

for i=1:9;
    ampsNear = nanmax(PeakTypeBurst{i},[],1);
    ampsFar = nanmax(PeakTypeOther{i},[],1);
    figure;
    hold on;
    [CdfY,CdfX,flo,fup] = ecdf(ampsNear','Function','cdf');
    CdfX=CdfX(2:end)';
    CdfY=CdfY(2:end)';
    flo=flo(2:end)';
    fup=fup(2:end)';
    X = [CdfX,fliplr(CdfX)];
    Y = [flo(1:end-1),CdfY(end),CdfY(end), fliplr(fup(1:end-1))];
    plot(CdfX,CdfY,'ok','MarkerSize',4,'MarkerFaceColor','k');
    patch(X,Y,'b','FaceAlpha',0.2,'EdgeColor','b');
    
    [CdfY,CdfX,flo,fup] = ecdf(ampsFar','Function','cdf');
    CdfX=CdfX(2:end)';
    CdfY=CdfY(2:end)';
    flo=flo(2:end)';
    fup=fup(2:end)';
    X = [CdfX,fliplr(CdfX)];
    Y = [flo(1:end-1),CdfY(end),CdfY(end), fliplr(fup(1:end-1))];
    plot(CdfX,CdfY,'sk','MarkerSize',4,'MarkerFaceColor','k');
    patch(X,Y,'r','FaceAlpha',0.2,'EdgeColor','r');
    legend({'','Near','','Far'});
    
    
end
