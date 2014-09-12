function [PeakTypeBurst,PeakTypeOther]=ClassifyAstroPeaks(t,ic,traces,time,bs,be)
% Function which finds the peaks in the astro trace which are associated with bursts. 
[PairWisePSTH,PairWiseLags]=CalcPSTHastrotriggers(t,ic,traces,time,bs,be);

%in progress
end