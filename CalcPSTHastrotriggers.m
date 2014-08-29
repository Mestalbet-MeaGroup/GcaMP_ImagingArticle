function [PairWisePSTH,PairWiseLags]=CalcPSTHastrotriggers(t,ic,traces,time,bursts,burste)

pre= 5000; %ms
post = pre;
bin = 100; %ms

% Create lags
[t,ic] = RemoveHAneurons(t,ic,bursts,burste);
temp  = mpsth(t(ic(3,1):ic(4,1))./12000,time(floor(end/2)),'tb',1,'fr',0,'pre',pre,'post',post,'binsz',bin);
PairWiseLags = temp(:,1);

% Detect amplitude rises
[StartAmpRise,~,~]=CalcPeakStartEnd(traces);
c=1;
% calculate pairwaise psth
tic
for i=1:size(traces,2)
    for j=1:size(ic,2)
        ps = mpsth(t(ic(3,j):ic(4,j))./12000,time(StartAmpRise{i}),'tb',0,'fr',0,'pre',pre,'post',post,'binsz',bin);
        PairWisePSTH(c,:) = ps;
        c=c+1;
    end
end
toc

tic
combs = allcombs(1:size(traces,2),1:size(ic,2));
i = combs(:,1);
j = combs(:,2);
starts = ic(3,j);
ends = ic(4,j);
fun = @(x,s,e) x(s:e);
tt = arrayfun(fun, t,starts,ends,'UniformOutput',false);
fun = @(x,s,i) x(s{i});
ttime = arrayfun(fun, time,StartAmpRise,i,'UniformOutput',false);

parfor c=1:size(combs,1)
    ps = mpsth(tt{c}./12000,time(StartAmpRise{i(c)}),'tb',0,'fr',0,'pre',pre,'post',post,'binsz',bin);
    PairWisePSTH(c,:) = ps;
end
toc

end
