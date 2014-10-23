function [b2pi,p2bi,p2b_cv,b2p_cv,peakIndex]=FindPeaksNearestBursts(traces,time,BurstOnsets,t);
% Function which take traces, times and burst starts and (1) Finds the
% peaks within each trace and (2) Calculates the minimum distance between
% each peak and its nearest burst start.

[~,peakIndex,~]=CalcPeakStartEnd(traces);
b = arrayfun(@(x) max(t(find(t>=x,10,'First'))),BurstOnsets)./12000; %convert to seconds and look at 10th spike as burst start (avoid premature burst start detection)

for i=1:size(traces,2) %ROI
    if ~isempty(peakIndex{i})
        a = time(peakIndex{i}); %time of peak in seconds
        [c,ai] = sort(b);
        [~,ic] = histc(a,[-inf,(c(1:end-1)+c(2:end))/2,inf]);
        p2bi{i} = ai(ic);
        p2b_cv{i} =a-b(ai(ic));
        
        [c,ai] = sort(a);
        [~,ic] = histc(b,[-inf,(c(1:end-1)+c(2:end))/2,inf]);
        b2pi{i} = ai(ic);
        b2p_cv{i} =b-a(ai(ic));
    end
end
% for i=1:size(peakIndex,2)
%     a = time(peakIndex{i});
%     for j=1:numel(b)
%         for k=1:numel(a)
%             dist(j,k)=a(k)-b(j);
%         end
%     end
% end
end