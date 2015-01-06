function [numBI,numNoB,numSB,PeaksWithBursts,BurstsWithPeaks,PeaksNoBursts,BurstsNoPeaks,PeaksWithSBs,SBsWithPeaks,peakIndex]=CalcBurstPeakOccurenceRate(DataSet,k);


[~,peakIndex,~]=CalcPeakStartEnd(DataSet{k}.dfTraces);
peaks = DataSet{k}.dfTime(cell2mat(peakIndex'));
distances = bsxfun(@minus,peaks',DataSet{k}.bs./12000);
closeby = zeros(size(distances));
closeby = abs(distances)<1;
numbs = sum(DataSet{k}.bs<=max(peaks)*12000);
validBS = find( (DataSet{k}.bs<=max(peaks)*12000) & (DataSet{k}.bs>=min(peaks)*12000) );
% closeby = closeby(:,validBS);
numbs=numel(validBS);
numBI = sum(sum(closeby,1)>=1)/numbs;

%---Reference Peaks to ROI---%
peaksperROI = cellfun(@(x) 1:numel(x),peakIndex,'UniformOutput',false);
add = cumsum(cellfun(@(x) numel(x),peakIndex));
add = num2cell(add);
peaksperROI(2:end) = cellfun(@(x,y) x+y,peaksperROI(2:end),add(1:end-1),'UniformOutput',false);
%---------------------------%
BurstsWithPeaks=cell(1,size(closeby,2));
for i=1:size(closeby,2)
    for j=1:size(peaksperROI,2)
        if sum(closeby(peaksperROI{j},i))>=1
            BurstsWithPeaks{i} = [BurstsWithPeaks{i};j];
        end
    end
end

for i=1:size(peaksperROI,2)
%     PeaksWithBursts{i} = find(sum(closeby(peaksperROI{i},:),1)>=1);
      [PeaksWithBursts{i},~] = find(closeby(peaksperROI{i},:)==1);
end

%%
if ~isempty(DataSet{k}.sbs)
    distances = bsxfun(@minus,peaks',DataSet{k}.sbs./12000);
    closeby = zeros(size(distances));
    closeby = abs(distances)<mean(DataSet{k}.sbe-DataSet{k}.sbs);
    numSB = sum(sum(closeby,1)>=1)/numel(DataSet{k}.sbs);
    distances = bsxfun(@minus,peaks',[DataSet{k}.bs./12000,DataSet{k}.sbs./12000]);
%     distances = bsxfun(@minus,peaks',DataSet{k}.bs./12000);
    
else
    numSB=nan;
end

SBsWithPeaks=cell(1,size(closeby,2));
for i=1:size(closeby,2)
    for j=1:size(peaksperROI,2)
        if sum(closeby(peaksperROI{j},i))>=1
            SBsWithPeaks{i} = [SBsWithPeaks{i};j];
        end
    end
end

for i=1:size(peaksperROI,2)
    PeaksWithSBs{i} = find(sum(closeby(peaksperROI{i},:),1)>=1);
end

%%


far = zeros(size(distances));
far = abs(distances)<=3;

numNoB = sum(sum(far,2)==0)/numel(peaks);

BurstsNoPeaks=cell(1,size(far,2));
for i=1:size(far,2)
    for j=1:size(peaksperROI,2)
        if sum(far(peaksperROI{j},i))==0
            BurstsNoPeaks{i} = [BurstsNoPeaks{i};j];
        end
    end
end


for i=1:size(peaksperROI,2)
    temp  = find(sum(far(peaksperROI{i},:),2)==0);
    if ~isempty(DataSet{k}.sbs)
        for j=1:numel(DataSet{k}.sbs)
            locs = find((DataSet{k}.dfTime(peakIndex{i}(temp))*12000)>=DataSet{k}.sbs(j) & (DataSet{k}.dfTime(peakIndex{i}(temp))*12000)<=DataSet{k}.sbe(j));
            temp(locs)=[];
        end
    end
    PeaksNoBursts{i} = temp;
end

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