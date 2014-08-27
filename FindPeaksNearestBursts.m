function [peakIndex,OnIndex,OffIndex,ClosestValues]=FindPeaksNearestBursts(traces,time,BurstOnsets);

f = savgol(35,5,0);
tr=filtfilt(f,1,traces);
for j=1:size(tr,2)
    peakIndex{j} = peakfinder(tr(:,j),(max(tr(:,j))-min(tr(:,j)))/4);
    for k=1:size(peakIndex{j},1)
        bs(k)  = RollDownBack(tr(:,j),peakIndex{j}(k),mean(tr(:,j))*1.05);
        betry = find(tr(peakIndex{j}(k):end,1)==0,1,'First')+peakIndex{j}(k);
        if isempty(betry)
            be(k)  = RollDownFront(tr(:,j),peakIndex{j}(k),mean(tr(:,j))*1.05);
        else
            be(k)=betry;
        end
    end
    sames = sort([find(diff(bs)==0),find(diff(bs)==0)+1]);
    [~,loc]=max(tr(sames,1));
    remove = setdiff(sames,sames(loc));
    peakIndex{j}(remove)=[];
    be(remove)=[];
    bs(remove)=[];
    OnIndex{j}=bs;
    OffIndex{j}=be;
    clear bs; clear be;
    a = time(OnIndex{j});
    b = BurstOnsets;
    m = numel(a); n = numel(b);
    [c,p] = sort([a,b]);
    q = 1:m+n; q(p) = q;
    t = cumsum(p>m);
    r = 1:n; r(t(q(m+1:m+n))) = r;
    s = t(q(1:m));
    id = r(max(s,1));
    iu = r(min(s+1,n));
    [d,it] = min([abs(a-b(id));abs(b(iu)-a)]);
    if var(d)==0
        ClosestValues(j)=nan;
    else
%         ClosestValues(j)=mean(d)/var(d);
        ClosestValues(j)=1/var(d);
    end
end


end