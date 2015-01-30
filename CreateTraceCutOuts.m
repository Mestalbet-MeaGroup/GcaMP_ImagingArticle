function cutout = CreateTraceCutOuts(trace,center,before,after)

start = center-before;
stop = center+after;
if start<1
    temp = padarray(trace,[abs(start),0],'replicate','pre');
    cutout = temp(center+abs(start)-before+1:center+abs(start)+after+1);
else
    if stop>numel(trace)
        offset = stop - numel(trace); 
        temp = padarray(trace,[offset,0],'replicate','post');
        cutout = temp(start:center+after);
    else
        cutout = trace(start:stop);
    end
end
end