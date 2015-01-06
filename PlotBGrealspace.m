function PlotBGrealspace(mask,connmat,emptynodes,noderank);
% Function which accepts the mask, connection matrix and the nodes to remove and produces a plot of the connections between astrocytes
% in real space.
opengl software
connmat=connmat+connmat';
connmat(emptynodes,:)=[];
connmat(:,emptynodes)=[];

% centrality = pagerank(connmat>0);
% [~,~,noderank] = unique(centrality,'stable');
nodecolor = parula(max(noderank));

pixels = regionprops(~mask,'PixelIdxList');
pixels(emptynodes)=[];

image = zeros(size(mask));
for i=1:numel(noderank)
    image(pixels(i).PixelIdxList)=noderank(i);
end

cdatao = real2rgb(image,[0,0,0;nodecolor],[0,max(noderank)]);
xo = repmat(1:size(mask,2),size(mask,1),1).*0.6472;
yo = repmat([1:size(mask,1)]',1,size(mask,2)).*0.6472;
zo = zeros(size(mask));
s = surface(xo,yo,zo,cdatao,'EdgeColor','none','FaceColor','texturemap','CDataMapping','direct');
axis('tight');
set(gca,'FontSize',18,'TickDir','out','PlotBoxAspectRatio',[size(mask,2),size(mask,1),1]);

centers  = cell2mat(arrayfun(@(x) x.Centroid, regionprops(~mask,'Centroid'),'uniformoutput',0));
centers(emptynodes,:)=[];
centers = centers.*0.6472;
[a1,a2]=find(connmat>0);
xposa1 = (centers(a1,1))';
yposa1 = (centers(a1,2))';
xposa2 = (centers(a2,1))';
yposa2 = (centers(a2,2))';
edgeval = arrayfun(@(x,y) connmat(x,y),a1,a2);
[~,~,edgerank] = unique(round(edgeval*1000)/1000,'stable');
edgecolor = parula(max(edgerank));

num2disp = 200;
randedges = randperm(numel(a1));
randedges = randedges(1:num2disp);
for i=1:num2disp
    l = line([xposa1(randedges(i)),xposa2(randedges(i))],[yposa1(randedges(i)),yposa2(randedges(i))],'linestyle','-','Color',edgecolor(edgerank(randedges(i)),:),'linewidth',1);
end
set(gca,'FontSize',9);

a = [a1,a2];
x = [xposa1,xposa2];
y = [yposa1,yposa2];
[~,ib,~]=unique(a,'stable');
a=a(ib);
x=x(ib);
y=y(ib);
a1labels =cellfun(@(x) num2str(x), num2cell(a),'UniformOutput',false);
t = text(x,y,a1labels,'Clipping', 'on','hittest','off','color','r','FontSize',18,'FontWeight','bold');
% arrayfun(@(x) uistack(x,'top'),t);
c1 = colorbar('Location','eastoutside');
% ctick = 0:(1/numel(avals))/2:1;
% ctick=ctick(2:2:end-1);
set(c1,'Limit',[min(edgerank),max(edgerank)],'Ticks',[min(edgerank),max(edgerank)],'TickLabels',{'lowest','highest'});
ylabel(c1,'rank')
set(gca,'FontSize',9,'TickDir','out');
end