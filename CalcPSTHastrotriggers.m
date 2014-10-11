function [PairWisePSTH,PairWiseLags]=CalcPSTHastrotriggers(t,ic,traces,time,bursts,burste)

pre= 5000; %ms
post = pre;
bin = 100; %ms

% Create lags
[t,ic] = RemoveHAneurons(t,ic,bursts,burste);
[t,ic] = CutSortChannel2(t,ic,min(time)*12000,max(t));
temp  = mpsth(t(ic(3,1):ic(4,1))./12000,time(floor(end/2)),'tb',1,'fr',0,'pre',pre,'post',post,'binsz',bin);
PairWiseLags = temp(:,1);

% Detect amplitude rises
[~,StartAmpRise,~]=CalcPeakStartEnd(traces); %Peaks
% [StartAmpRise,~,~]=CalcPeakStartEnd(traces); %Peak starts

combs = allcomb(1:size(traces,2),1:size(ic,2));
assignin('base','combs',combs);
i = combs(:,2);
j = combs(:,1);
starts = ic(3,i);
ends = ic(4,i);
StartAmpRise=StartAmpRise(j);
fun = @(s,e) t(s:e);
tt = arrayfun(fun, starts,ends,'UniformOutput',false);

fun = @(x) time(x);
ttime = cellfun(fun, StartAmpRise,'UniformOutput',false);

parfor c=1:size(combs,1)
    ps = mpsth(tt{c}./12000,ttime{c},'tb',0,'fr',1,'pre',pre,'post',post,'binsz',bin)./numel(ttime{c});
    PairWisePSTH(c,:) = ps;
end
% calculate pairwaise psth serially
% for i=1:size(traces,2)
%     for j=1:size(ic,2)
%         ps = mpsth(t(ic(3,j):ic(4,j))./12000,time(StartAmpRise{i}),'tb',0,'fr',0,'pre',pre,'post',post,'binsz',bin);
%         PairWisePSTH(c,:) = ps;
%         c=c+1;
%     end
% end
end
