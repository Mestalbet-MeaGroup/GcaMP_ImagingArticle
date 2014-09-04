function s = Plot2DHist(x,y,binx,biny,xlimit)
% Function which takes two dependent variables (like distance and
% max correlation or lag and max correlation) and plots a 2d histogram
% where color scales with number of occurence
[hist,bins]  = hist3([x,y],[binx,biny]);
xo = repmat(bins{1}',1,size(hist,2));
yo = repmat(bins{2},size(hist,1),1);
zo = zeros(size(hist));
[rank,~,~] = unique(hist(:));
colors = [1,1,1;jet(numel(rank))];
cdatao = real2rgb(hist,colors,[0,max(rank)]);
s = surface(xo,yo,zo,cdatao,'EdgeColor','none','FaceColor','texturemap','CDataMapping','direct');
% set(gca,'XLim',[0,3100],'Ylim',[0,1]);
set(gca,'XLim',[0,xlimit],'Ylim',[0,1]);
colormap(colors);
colorbar;

%---Center of Mass---%
[rc,cc] = ndgrid(1:size(hist,1),1:size(hist,2));
Mt = sum(hist(:));
c1 = sum(hist(:) .* rc(:)) / Mt;
c2 = sum(hist(:) .* cc(:)) / Mt;
x=interp1(1:size(bins{1},2),bins{1},c1);
y=interp1(1:size(bins{2},2),bins{2},c2);
hold on;
plot(x,y,'+k','MarkerSize',5,'LineWidth',10);
end