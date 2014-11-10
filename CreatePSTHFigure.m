function CreatePSTHFigure(DataSet,lag,psth,i)
subplot = @(m,n,p) subtightplot (m, n, p, [0.025 0.025], [0.1 0.05], [0.05 0.1]);

figure('Color','White');
s1 = subplot(5,1,1:3);
cult = DataSet{i,1}.culture;
ch = DataSet{i,1}.channel;
title(['Astrocyte Triggered, ','Culture ',num2str(cult),' Channel: ',num2str(ch)]);
%---Order by Offsets---%
% order = Sortpsth{i}pairs(atrigpsth{i}{i});
%---No Ordering---%
order = 1:size(psth{i},1);
%----------------------%
temp = zscore(psth{i}(order,:)')';
temp(temp<0)=0;
temp = temp(2500:end,:); 
imagesc(lag{i},1:size(temp,1),temp);
xlim([-4900,4900]);
ylabel('Pairs');
set(gca,'TickDir','out','XTick',[],'TickLength',[0.0025,0.0025]);
hb = colorbar('location','eastoutside');
set(hb,'TickDir','out');
s1Pos = get(s1,'position');
set(gca,'FontSize',18,'FontName','Times-Roman','YDir','normal');

s2 = subplot(5,1,4);
errorbar(lag{i},nanmean(zscore(psth{i},0,2),1),nanstd(zscore(psth{i},0,2),[],1),'.-')
axis tight;
xlim([-4900,4900]);
set(gca,'TickDir','out','XTick',[]);
s2Pos = get(s2,'position');
s1Pos(3) = s2Pos(3);
set(s1,'position',s1Pos);
box off;
ylabel('<Z-Score>');
set(gca,'FontSize',18);

s3 = subplot(5,1,5);
atrigmean =[];
hold all
for i=1:9
    plot(lag{i},nanmean(zscore(psth{i},0,2),1));
    atrigmean = [atrigmean;nanmean(zscore(psth{i},0,2))];
end
plot(lag{i},nanmean(atrigmean,1),'--k','LineWidth',3);
axis tight;
xlim([-4900,4900]);
set(gca,'TickDir','out');
legs = cellfun(@(x) ['Network ',num2str(x)],num2cell(1:size(psth,2)),'UniformOutput',0);
legs(end+1)={'Mean'};
legend(legs,'location','eastoutside');
s3Pos = get(s3,'position');
s3Pos(3) = s2Pos(3);
set(s3,'position',s3Pos);
xlabel('time [ms]');
ylabel('<Z-Score>');
set(gca,'FontSize',18,'FontName','Times-Roman');

end
