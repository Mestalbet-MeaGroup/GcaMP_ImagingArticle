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

% calculate pairwaise psth
for i=1:size(traces,2)
    for j=1:size(ic,2)
        ps = mpsth(t(ic(3,j):ic(4,j))./12000,time(StartAmpRise{i}),'tb',0,'fr',0,'pre',pre,'post',post,'binsz',bin);
        PairWisePSTH(i,j,:) = ps;   
    end
end