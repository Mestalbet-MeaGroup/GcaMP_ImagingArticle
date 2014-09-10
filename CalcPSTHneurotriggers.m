function [PairWisePSTH,PairWiseLags]=CalcPSTHneurotriggers(t,ic,traces,time,bursts,burste)

pre= 5000; %ms
post = pre;
bin = 100; %ms

% Create lags
% [t,ic] = RemoveHAneurons(t,ic,bursts,burste);
temp  = mpsth(t(ic(3,1):ic(4,1))./12000,time(floor(end/2)),'tb',1,'fr',0,'pre',pre,'post',post,'binsz',bin);
PairWiseLags = temp(:,1);

% Detect amplitude rises
% [~,pp,~]=CalcPeakStartEnd(traces); %peaks
[pp,~,~]=CalcPeakStartEnd(traces); %peak starts
%-----%

fun = @(x) time(x);
tt = cellfun(fun, pp,'UniformOutput',false);
PairWisePSTH=zeros(size(pp,2),size(PairWiseLags,1));
parfor c=1:size(pp,2)
    ps = mpsth(tt{c},bursts./12000,'tb',0,'fr',0,'pre',pre,'post',post,'binsz',bin)./numel(tt{c});
    PairWisePSTH(c,:) = ps;
end

end
