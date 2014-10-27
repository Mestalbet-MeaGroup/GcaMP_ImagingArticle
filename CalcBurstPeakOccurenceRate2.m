function [numBI,numNoB,numSB]=CalcBurstPeakOccurenceRate2(DataSet,k);


[~,peakIndex,~]=CalcPeakStartEnd(DataSet{k}.dfTraces);
peaks = DataSet{k}.dfTime(unique(cell2mat(peakIndex')));
distances = bsxfun(@minus,peaks',DataSet{k}.bs./12000);
closeby = zeros(size(distances));
closeby = abs(distances)<1;
numbs = sum(DataSet{k}.bs<=peaks(end)*12000);
numBI = sum(sum(closeby,1)>=1)/numbs;
%%
if ~isempty(DataSet{k}.sbs)
    distances = bsxfun(@minus,peaks',DataSet{k}.sbs./12000);
    closeby = zeros(size(distances));
    closeby = abs(distances)<3;
    numSB = sum(sum(closeby,1)>=1)/numel(DataSet{k}.sbs);
    distances = bsxfun(@minus,peaks',[DataSet{k}.bs./12000,DataSet{k}.sbs./12000]);
else
    numSB=nan;
end
%%
far = zeros(size(distances));
far = abs(distances)<=3;
numNoB = sum(sum(far,2)==0)/numel(peaks);

%% Plot
% hold on;
% h= bar(bars,'stacked');
% set(h(1),'facecolor',[84,95,255]./255);
% set(h(2),'facecolor',[157,118,208]./255);
% 
% barsx = get(h,'XData');
% barsy = get(h,'YData');
% 
% h=errorbar(barsx{1},barsy{1},ebars(:,1),'.r');
% % h=errorbar(barsx{2},barsy{2},ebars(:,2),'.');
% axis tight;
% set(gca,'XTickLabel',{'Bursts With Peaks','','Superbursts With Peaks','','Peaks Without Bursts'});
% set(gca,'FontSize',18,'TickDir','Out');
% ylabel('Occurence');
% ylim([0,1.1]);




end