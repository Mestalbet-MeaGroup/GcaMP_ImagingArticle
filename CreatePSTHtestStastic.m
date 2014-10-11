function [Rpsth,lags] = CreatePSTHtestStastic(t,ic,traces,time,bursts,burste);

pre= 5000; %ms
post = pre;
bin = 100; %ms
% Create lags
[t,ic] = RemoveHAneurons(t,ic,bursts,burste);
[t,ic] = CutSortChannel2(t,ic,min(time)*12000,max(t));
temp  = mpsth(t(ic(3,1):ic(4,1))./12000,time(floor(end/2)),'tb',1,'fr',0,'pre',pre,'post',post,'binsz',bin);
lags = temp(:,1);

% Detect amplitude rises
[~,StartAmpRise,~]=CalcPeakStartEnd(traces); %Peaks

combs = allcomb(1:size(traces,2),1:size(ic,2));
i = combs(:,2);
j = combs(:,1);
starts = ic(3,i);
ends = ic(4,i);
StartAmpRise=StartAmpRise(j);
fun = @(x) time(x);
ttime = cellfun(fun, StartAmpRise,'UniformOutput',false);

for k=1:10
   psth = CalcRandPSTH(t,ic,starts,ends,ttime,size(combs,1),numel(lags));
   Rpsth(k,:)=nanmean(psth,1);
end


end