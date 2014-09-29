function PlotCentralityRealSpace(e1,pr,a1,a_e2,avals,pos,viewfield)

% f = figure('Visible','off');
hold on;

%---Load Mea Image Data and Map---%
load('MeaMapPlot.mat','MeaImage','MeaMap');
MeaMap = rot90(MeaMap,-1);
overlay = otsu(imresize(MeaImage,1.0441),45);
overlay(overlay>1)=0;
overlay=bwareaopen(overlay, 500);
overlay = (overlay>=1);
%--Erase manually two spots---%
overlay(1750:1830,2705:2770)=0;
overlay(1600:1720,2745:2830)=0;
%-----------------------------%
centers = imfindcircles(overlay,[12,35],'ObjectPolarity','bright','method','TwoStage','Sensitivity',0.70);
[xpos,a]=sort(centers(:,2));
ypos=centers(a,1);

%--Find Electrode Positions---%
mask = otsu(MeaImage,45);
mask(mask>1)=0;
centers = imfindcircles(mask,[12,20],'ObjectPolarity','bright','method','TwoStage','Sensitivity',0.90);
[xpos,a]=sort(centers(:,2));
ypos=centers(a,1);
%---Add Ground Electrodes to x and y pos vectors---%
xpos = [nan;xpos(1:14);nan;xpos(15:238);nan;xpos(239:end);nan];
ypos = [nan;ypos(1:14);nan;ypos(15:238);nan;ypos(239:end);nan];

%---Reshape xpos and ypos into 16x16 grid---%
xpos=reshape(xpos,16,16);
ypos=reshape(ypos,16,16);
%--Plot Mask---%
[sitex,sitey]=find(MeaMap==pos);
viewfield = imresize(viewfield,[357.6,501.2]);
cdatao = real2rgb(viewfield,[255,0,0;152,152,152]./255,[0,1]);
xo = repmat(1:size(viewfield,2),size(viewfield,1),1)+xpos(sitex,sitey)-2*size(viewfield,1)/3;
yo = repmat([1:size(viewfield,1)]',1,size(viewfield,2))+ypos(sitex,sitey)-size(viewfield,2)/3-12.5;
zo = zeros(size(viewfield));
hf=figure('color','white');
whitebg('black')
hold on;
s = surface(xo,yo,zo,cdatao,'EdgeColor','none','FaceColor','texturemap','CDataMapping','direct'); axis('tight');%axis('square');
xvec = 0:480.3997:2.8824e+03;
set(gca,'XTick',xvec,'XTickLabel',xvec.*1.0408); % 200/sqrt((xpos(1,3)-xpos(1,2))^2+(ypos(1,3)-ypos(1,2))^2)
uistack(s,'top');
alpha(s,'texture');


[uPr,~,ranks]=unique(pr);
colors = jet(numel(uPr));

msize = round(((ranks-min(ranks))/(max(ranks)-min(ranks))).*30)+10;
plot(centers(:,2),centers(:,1), 'o',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1 1 1]);
for i=1:size(pr)
    [x1,y1] = find(MeaMap==e1(i));
    plot(xpos(x1,y1),ypos(x1,y1),'.','color',colors(ranks(i),:),'MarkerSize',msize(i));
end
xlim([-100,3100]);
ylim([-100,3100]);

centers  = cell2mat(arrayfun(@(x) x.Centroid, regionprops(~viewfield,'Centroid'),'uniformoutput',0));
xposa = (centers(a1,1)+xo(1,1))';
yposa = (centers(a1,2)+yo(1,1))';
[allvals,rank]=sort(avals);
sources = a1;
colors = flip(jet(numel(sources)),1);
[~,~,allv_index] = unique(allvals);
for i=1:numel(a_e2)
    x1=xposa(i);
    y1=yposa(i);
    [ax2,ay2] = find(MeaMap==a_e2(i));
    loc = find(sources==(a1(i)),1,'First');
    score=colors(rank(loc),:);
    sources(loc)=nan; %in case two connections have the same source (non-unique elements in a1)
    al = line([x1,xpos(ax2,ay2)],[y1,ypos(ax2,ay2)],'color',score,'linewidth',0.5);
    uistack(al,'top');
end

end