function [WhichTraces,start,stop] = CalcBestTrio(pi,traces,time,secs,fs);

window=ceil(secs*fs);
traces = zscore(traces);
numbers=[];
pkidx=[];
for i=1:size(pi,2)
    if ~isempty(pi{i})
        numbers=[numbers,ones(1,numel(pi{i})).*i];
        pkidx = [pkidx,pi{i}];
    end
end
peaks = [numbers;pkidx]';
for i=1:size(peaks,1)
    peaks(i,3)=traces(pkidx(i),numbers(i));
end
peaks = sortrows(peaks,2);

combs = VChooseK(1:size(peaks,1),3);
c1 = combs(:,1);
c2 = combs(:,2);
c3 = combs(:,3);
p1a = peaks(c1,1);
p2a = peaks(c2,1);
p3a = peaks(c3,1);
p1b = peaks(c1,2);
p2b = peaks(c2,2);
p3b = peaks(c3,2);
p1c = peaks(c1,3);
p2c = peaks(c2,3);
p3c = peaks(c3,3);
c1 =zeros(numel(combs(:,1)),1);
c2 =zeros(numel(combs(:,1)),1);
c3 =zeros(numel(combs(:,1)),1);
parfor i=1:numel(c1)
    c1(i) = numel(unique([p1a(i),p2a(i),p3a(i)]));
    temp1  = sort([p1b(i),p2b(i),p3b(i)]);
    c2(i) = temp1(end) - temp1(1);
    c3(i) = median([p1c(i),p2c(i),p3c(i)])/(max([p1c(i),p2c(i),p3c(i)])-min([p1c(i),p2c(i),p3c(i)]));
end

get3 = find(c1==3);
getclose = find(c2<ceil(1*fs));
gethigh = find(c3>5);
runIntersect = mintersect(get3,getclose,gethigh);
if isempty(runIntersect)
    get3 = find(c1==3);
    getclose = find(c2<ceil(2*fs));
    gethigh = find(c3>5);
    runIntersect = mintersect(get3,getclose,gethigh);
end
if isempty(runIntersect)
    get3 = find(c1==3);
    getclose = find(c2<ceil(2*fs));
    gethigh = find(c3>5);
    runIntersect = mintersect(get3,getclose,gethigh);
end
if isempty(runIntersect)
    get3 = find(c1==3);
    getclose = find(c2<ceil(15*fs));
    gethigh = find(c3>3);
    runIntersect = mintersect(get3,getclose,gethigh);
end
candidate = find(max(c3(runIntersect))==c3);
WhichTraces= peaks(combs(candidate,:),1);
extra=0;
start = floor(median(peaks(combs(candidate,:),2)))-ceil(window/2);
if start<=0
    extra = abs(start);
    start=1;
end
start = time(start)*12000;
stop = time(floor(median(peaks(combs(candidate,:),2)))+ceil(window/2)+extra)*12000;
end