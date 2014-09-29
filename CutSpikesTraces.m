function [t,ic,traces,time] = CutSpikesTraces(t,ic,bs,be,traces,time,start,stop);
%Function which takes a start time in seconds and an end time in seconds
%and cuts the data within this window.
[t,ic,bs,be,~] = SortChannelsByFR2(t,ic,bs,be,be-bs);
[t,ic] = OnlyNonHA(t,ic,bs,be);
[t,ic]=CutSortChannel2(t,ic,start*12000,stop*12000);
index = (time>=start)&(time<=stop);
time = time(index);
traces = traces(index,:);
time = time - min(time);
end