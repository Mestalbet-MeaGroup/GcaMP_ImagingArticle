function [CutoutNear, CutoutFar, NearROIs, FarROIs, ampvDist] = CalcCutOuts(DataSet,k,counter,cf)
    
%     [~,peakIndex,~]=CalcPeakStartEnd(DataSet{k}.dfTraces);
[~,~,~,PeaksWithBursts,~,PeaksNoBursts,~,~,~,peakIndex]=CalcBurstPeakOccurenceRate(DataSet,k);
for i=1:size(peakIndex,2)
peaks = DataSet{k}.dfTime(peakIndex{i});
distances = bsxfun(@minus,peaks',DataSet{k}.bs./12000);
ampvDist{i} = min(abs(distances),[],2);
tr  = zscore(DataSet{k}.dfTraces(:,i));
ampvDist{i}(:,2)=tr(peakIndex{i});

%         closeby = zeros(size(distances));
%         closeby = abs(distances)<1;
%         far = zeros(size(distances));
%         far = abs(distances)<=3;
%
%         [pc{i},~]=find(closeby==1);
%         pf{i}=find(sum(far==0,2)==size(distances,2));

if ~isempty(PeaksWithBursts{i})
    PeaksWithBursts{i}=unique(PeaksWithBursts{i});
    for j=1:numel(PeaksWithBursts{i})
        CutoutNear(counter,:) = CreateTraceCutOuts(tr,peakIndex{i}(PeaksWithBursts{i}(j)),50,50);
        %                 CutoutNear(counter,:) = CreateTraceCutOuts(DataSet{k}.dfTraces(:,i),peakIndex{i}(pc{i}(j)),50,50);
        NearROIs(counter,:) = [k,PeaksWithBursts{i}(j)];
        counter=counter+1;
    end
end

if ~isempty(PeaksNoBursts{i})
    PeaksNoBursts{i}=unique(PeaksNoBursts{i});
    for j=1:numel(PeaksNoBursts{i})
        CutoutFar(cf,:) = CreateTraceCutOuts(tr,peakIndex{i}(PeaksNoBursts{i}(j)),50,50);
        %                 CutoutFar(cf,:) = CreateTraceCutOuts(DataSet{k}.dfTraces(:,i),peakIndex{i}(pf{i}(j)),50,50);
        FarROIs(cf,:) = [k,PeaksNoBursts{i}(j)];
        cf=cf+1;
    end
end
end

end