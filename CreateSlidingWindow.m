function [newx,newy]=CreateSlidingWindow(xvals,yvals,nwind,noverlap)
[xvals,index]=sort(xvals);
yvals=yvals(index);

nx = length(xvals);                            %# length of sequence
ncol = fix((nx-noverlap)/(nwind-noverlap));    %# number of sliding windows
colindex = 1 + (0:(ncol-1))*(nwind-noverlap);  %# starting index of each
idx = bsxfun(@plus, (1:nwind)', colindex)-1;   %# each column is the indexes of one window. Each row is a sliding step forward. 
newx = nanmean(xvals(idx),1);                  %# apply the indices on the sequence and take the average
% newy = nanmean(yvals(idx),1);                  %# apply the indices on the sequence and take the average     
%---For All Values---%
% for i=1:size(idx,2)
%     newy{i}=yvals(idx(:,i));
% end
newy=yvals(idx);
end