function [ps,pp,pe]=CalcPeakStartEnd(traces)
set(0,'RecursionLimit',1000);
f = savgol(30,5,0);
tr=filtfilt(f,1,traces);
% tr=traces;
for j=1:size(tr,2)
    pp{j} = peakfinder(tr(:,j),(max(tr(:,j))-min(tr(:,j)))/2);
    for k=1:size(pp{j},1)
        ps{j}(k)  = RollDownBack(tr(:,j),pp{j}(k),mean(tr(:,j))*1.05);
        betry = find(tr(ps{j}(k):end,1)==0,1,'First')+ps{j}(k);
        if isempty(betry)
            pe{j}(k)  = RollDownFront(tr(:,j),pp{j}(k),mean(tr(:,j))*1.05);
        else
            pe{j}(k)=betry;
        end
    end
    sames = sort([find(diff(ps{j})==0),find(diff(ps{j})==0)+1]);
    [~,loc]=max(tr(sames,1));
    remove = setdiff(sames,sames(loc));
    pp{j}(remove)=[];
    ps{j}(remove)=[];
    pe{j}(remove)=[];
end
end