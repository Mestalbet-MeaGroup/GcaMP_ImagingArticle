function s = Plot2DHist(x,y,binx,biny,xlimit,ylimit,factor,bins)
% Function which takes two dependent variables (like distance and
% max correlation or lag and max correlation) and plots a 2d histogram
% where color scales with number of occurence

interp=0;
%----Calculate Weighted Histogram---%
hist=zeros(binx,biny);
hist1=zeros(binx,biny);
hist2=zeros(binx,biny);
[~,idx]=histc(x,bins{1});
[~,idy]=histc(y,bins{2});
for i=1:size(idx,1)
    hist(idx(i)+1,idy(i)+1)=hist(idx(i)+1,idy(i)+1)+factor(i,idx(i)+1);
    hist1(idx(i)+1,idy(i)+1)=hist1(idx(i)+1,idy(i)+1)+1;
end

for k=1:binx
    f = trapz(bins{2},hist1(k,:));
    f(f==0)=inf;
    hist2(k,:) = hist(k,:)./f;
end

%---Compare Normalizations---%
figure;
subplot(1,3,1); imagesc(hist1'); colorbar; %set(gca,'PlotBoxAspectRatio',[binx,biny,1]);
subplot(1,3,2); imagesc(hist');  colorbar; %set(gca,'PlotBoxAspectRatio',[binx,biny,1]);
subplot(1,3,3); imagesc(hist2'); colorbar; %set(gca,'PlotBoxAspectRatio',[binx,biny,1]);

hist = hist2;
figure;
%---Calculate Raw Histogram---%
% [hist,bins]  = hist3([x,y],[binx,biny]);
%---Plot 2D Histogram---%
xo = repmat(bins{1}',1,size(hist,2));
yo = repmat(bins{2},size(hist,1),1);
zo = zeros(size(hist));
if interp == 1
    hist(hist==0)=nan;
    hist=inpaint_nans(hist,1);
    hist = round(hist.*100)./100;
    [rank,~,~] = unique(hist(:));
    colors = jet(numel(rank));
else
    [rank,~,~] = unique(hist(:));
    colors = [1,1,1;jet(numel(rank))];
end
cdatao = real2rgb(hist,colors,[0,max(rank)]);
hold on;
s = surface(xo,yo,zo,cdatao,'EdgeColor','none','FaceColor','texturemap','CDataMapping','direct');
% set(gca,'XLim',[0,3100],'Ylim',[0,1]);
% set(gca,'XLim',[0,xlimit],'Ylim',[0,1]);
axis tight;
ylims = get(gca,'YLim');
set(gca,'XLim',[0,xlimit+10],'YLim',[0,ylims(2)]);
% set(gca,'XLim',[0,xlimit],'Ylim',[0,ylimit]);


%---Center of Mass---%
[rc,cc] = ndgrid(1:size(hist,1),1:size(hist,2));
Mt = sum(hist(:));
c1 = sum(hist(:) .* rc(:)) / Mt;
c2 = sum(hist(:) .* cc(:)) / Mt;
x1=interp1(1:size(bins{1},2),bins{1},c1);
y1=interp1(1:size(bins{2},2),bins{2},c2);
plot(x1,y1,'+k','MarkerSize',5,'LineWidth',10);

%---Colorbar---%
ha  = axes('visible','off','CLim', [1, size(colors,1)]);
colormap(ha,colors);
cbh = colorbar('peer', ha, 'v',...
    'XTickLabel',{'0',num2str(nanmax(hist(:)))},...
    'XTick', [1,size(colors,1)],...
    'Position', [0.915077251802069,0.103321033210332,0.0138888888888888,0.81549815498155]);
end