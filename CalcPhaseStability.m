function ps=CalcPhaseStability(t,ic,varargin)
% function which takes the spike times (in any unit) and index channel, and calculates
% the normalized standard deviation of differences between closest spikes
% across all channels.
% Inputs:
% t, vector of spike times organized by channel
% ic, indexchannels (first row = channel id, second row = neuron id (relevant only in case of spike sorting), third
% row = index of t corresponding to first spike of channel, fourth row =
% index of t corresponding to last spike recorded on the channel)
% Optional Input:
% vector the same length as t, where each element is the channel the spike
% was recorded on (aka, Samora's method).
% Possible Usages:
% ps=CalcPhaseStability(t,[],ConvertIC2Samora(ic));
% or
% ps=CalcPhaseStability(t,ic);
% or
% ps=CalcPhaseStability(t,[],channelids);

if ~isempty(varargin)
    vec = varargin{1};
else
    vec = ConvertIC2Samora(ic);
end

list=unique(vec);
temp = VChooseKR(1:numel(list),2);
set1=list(temp(:,1));
set2=list(temp(:,2));
ps=zeros(numel(list),numel(list));
for i=1:numel(set1)
    diffs = FindMinDifferences(t(vec==set1(i)),t(vec==set2(i)));
%     [diffs,~,~] = deleteoutliers(diffs, 0.05);
%     maxval = abs(diffs-mean(diffs));
%     maxval=max(maxval);
%     if maxval > 0
%         ps(temp(i,1),temp(i,2))= 1 - std(diffs)/maxval;
%     else
%         ps(temp(i,1),temp(i,2))= 1 - std(diffs);
%     end
%     ps(temp(i,1),temp(i,2))= std(diffs)/mean(diffs);
    ps(temp(i,1),temp(i,2))=std(diffs);
end
ps=ps+ps'-eye(size(ps));
end

% must make this with a sliding window since a neuron can change who its
% locked to over the course of the recording. ps(:,:,time).