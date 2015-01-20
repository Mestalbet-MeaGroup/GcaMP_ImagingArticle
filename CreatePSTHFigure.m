function CreatePSTHFigure(DataSet,lag,psth,i,CutoutNear,ranksC)
subplot = @(m,n,p) subtightplot (m, n, p, [0.025 0.025], [0.1 0.05], [0.05 0.1]);

figure('Color','White');
minx=-3000;
maxx=2000;
s = subplot(6,1,1);
PlotPeaksPSTHplot(CutoutNear,ranksC,DataSet{i}.fs,minx/1000,maxx/1000);
s1 = subplot(6,1,2:4);
cult = DataSet{i,1}.culture;
ch = DataSet{i,1}.channel;
title(['Astrocyte Triggered, ','Culture ',num2str(cult),' Channel: ',num2str(ch)]);
%---Order by Offsets---%
order = SortPSTHpairs(zscore(psth{i},[],2));
%---No Ordering---%
% order = 1:size(psth{i},1);
%----------------------%
temp = zscore(psth{i}(order,:)')';
temp(temp<0)=0;
temp = temp(2500:end,:); 
imagesc(lag{i},1:size(temp,1),temp);
%xlim([-4900,4900]);
xlim([minx,maxx]);
ylabel('Pairs');
set(gca,'TickDir','out','XTick',[],'TickLength',[0.0025,0.0025]);
hb = colorbar('location','eastoutside');
set(hb,'TickDir','out');
s1Pos = get(s1,'position');
set(gca,'FontSize',18,'FontName','Times-Roman','YDir','normal');

s2 = subplot(6,1,5);
errorbar(lag{i},nanmean(zscore(psth{i},0,2),1),nanstd(zscore(psth{i},0,2),[],1),'.-')
axis tight;
% xlim([-4900,4900]);
xlim([minx,maxx]);
set(gca,'TickDir','out','XTick',[]);
s2Pos = get(s2,'position');
s1Pos(3) = s2Pos(3);
set(s1,'position',s1Pos);
box off;
ylabel('z-Score');
set(gca,'FontSize',18);

s3 = subplot(6,1,6);
atrigmean =[];
hold all
for i=1:9
    plot(lag{i},nanmean(zscore(psth{i},0,2),1));
    atrigmean = [atrigmean;nanmean(zscore(psth{i},0,2))];
end
plot(lag{i},nanmean(atrigmean,1),'--k','LineWidth',3);
axis tight;
%xlim([-4900,4900]);
xlim([minx,maxx]);
set(gca,'TickDir','out');
legs = cellfun(@(x) ['Network ',num2str(x)],num2cell(1:size(psth,2)),'UniformOutput',0);
legs(end+1)={'Mean'};
legend(legs,'location','eastoutside');
s3Pos = get(s3,'position');
s3Pos(3) = s2Pos(3);
set(s3,'position',s3Pos);
xlabel('time [ms]');
ylabel('z-Score');
set(gca,'FontSize',18,'FontName','Times-Roman');

pos(1,:) = get(s,'position');
for i=1:3
pos(i+1,:) = eval(['get(s',num2str(i),',''position'' )']);
end
h = axes('position',[pos(1,1),min(pos(:,2)),pos(1,3),sum(pos(:,4))+0.1208],'Color','none');
plot(h, zeros(1,2),[0,1],'Color',[204,0,102]./255,'LineStyle','--','LineWidth',5);
xlim(h,[minx,maxx]);
set(h,'Color','none');
axis(h,'off');
end

    
