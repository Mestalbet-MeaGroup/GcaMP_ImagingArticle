function PlotResultsOnMEAwithVF(e1,e2,value,varargin)
% Function which accepts a series of electrode sources and targets, and a
% value (same number of sources, targets and values) then plots a line
% connecting source to target, with a color determined by value. Requires
% the MEA map and MEA image mosaic files in MeaMapPlot.Mat.
% Optional Input:
% 1\ mask of ROIs (520x692)
% 2\ Viewfield location (electrod upon which viewfield is centered)
% 3\ Astrocytes ID (number in traces)
% 4\ Electrodes connected to astro's (must be same number as astrocytes (3))
% 5\ Magnitude of "connection" between astrocytes and electrodes
% Examples:
% To plot N2N (neuron to neuron):
% PlotResultsOnMEA(e1,e2,eval); where e1 is the source electrodes, e2 the targets and eval the value of the connection metric between them.
% To PlotA2N:
% PlotResultsOnMEA([],[],[],mask,ch,a1,a_e2,aval);
% To Plot both A2N and N2N:
% PlotResultsOnMEA(e1,e2,eval,mask,ch,a1,a_e2,aval);
% f = figure;
% Example call:
% [e1,e2,eval,a1,a_e2,aval,ch,mask]=ParseConnectionValuesWithVF(6,500);
% PlotResultsOnMEAwithVF(e1,e2,eval,mask,ch,a1,a_e2,aval);
%---Load Mea Image Data and Map---%
load('MeaMapPlot.mat','MeaImage','MeaMap');
MeaMap = rot90(MeaMap,-1);

%---Plot Mea Image---%

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
hold on
plot(centers(:,2),centers(:,1), 'o',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[1 1 1]);
xlim([-100,3100]);
ylim([-100,3100]);
set(gca,'XTick',[-100:250:3100],'XTickLabel',[-100:250:3100]+100);
set(gca,'YTick',[-100:250:3100],'YTickLabel',[-100:250:3100]+100);
set(gca,'Color','black','TickDir','out');
set(gca,'PlotBoxAspectRatio',[1,1,1]);
set(gca,'FontSize',18);

%--Plot Mask--%
if ~isempty(varargin)
    overlay = varargin{1};
    if size(varargin,2)>1
        pos = varargin{2};
        [sitex,sitey]=find(MeaMap==pos);
    else
        sitex=1;
        sitey=1;
    end
    mask = overlay;
    overlay = imresize(overlay,[357.6,501.2]);
%     cdatao = real2rgb(overlay,[255,0,0;152,152,152]./255,[0,1]);
    xo = repmat(1:size(overlay,2),size(overlay,1),1)+xpos(sitex,sitey)-2*size(overlay,1)/3;
    yo = repmat([1:size(overlay,1)]',1,size(overlay,2))+ypos(sitex,sitey)-size(overlay,2)/3-12.5;
    zo = zeros(size(overlay));
%     s = surface(xo,yo,zo,cdatao,'EdgeColor','none','FaceColor','texturemap','CDataMapping','direct'); %axis('tight');%axis('square');
%     xvec = 0:480.3997:2.8824e+03;
%     set(gca,'XTick',xvec,'XTickLabel',xvec.*1.0408); % 200/sqrt((xpos(1,3)-xpos(1,2))^2+(ypos(1,3)-ypos(1,2))^2)
%     uistack(s,'top');
%     alpha(s,'texture');
end

%--Find ROIs and Plot---%
if size(varargin,2)>2
    x1=[]; y1=[]; x2=[]; y2=[];ax2=[];ay2=[];
    a1=varargin{3};
    a_e2=varargin{4};
    avals=varargin{5};
    cM  = cell2mat(arrayfun(@(x) x.Centroid, regionprops(~mask,'Centroid'),'uniformoutput',0));
    centers=zeros(size(cM));
    [centers(:,1),centers(:,2)] = ConvertResizedCoordinates(mask,overlay,cM(:,1),cM(:,2));
%     centers  = cell2mat(arrayfun(@(x) x.Centroid, regionprops(~overlay,'Centroid'),'uniformoutput',0));
    xposa = (centers(a1,1)+xo(1,1))';
    yposa = (centers(a1,2)+yo(1,1))';
    [allvalsN,rank]=sort(value);
    [~,~,allv_index] = unique(allvalsN);
    sources = e1;
    sources=sources(rank);
   
    rank = rank(allv_index);
    sN=allvalsN;
    colorsN = cbrewer('seq','Purples',numel(sources));
    
    for i=1:numel(e1)
        [x1,y1] = find(MeaMap==e1(i));
        [x2,y2] = find(MeaMap==e2(i));
        loc = find(sources==e1(i),1,'First');
        score=colorsN(rank(loc),:);
        sources(loc)=nan;  %in case two connections have the same source (non-unique elements in e1)
%         l = patchline([xpos(x1,y1),xpos(x2,y2)],[ypos(x1,y1),ypos(x2,y2)],'linestyle','-','edgecolor',score,'linewidth',6,'edgealpha',0.8);
        l = line([xpos(x1,y1),xpos(x2,y2)],[ypos(x1,y1),ypos(x2,y2)],'color',score,'linewidth',1);
        uistack(l,'top');
    end
    [allvals,rank]=sort(avals);
    [~,~,allv_index] = unique(allvals);
    sA = allvals;
    sources = a1;
    sources=sources(rank);
    rank = rank(allv_index);
    colors = cbrewer('seq','YlOrRd',numel(sources));
    for i=1:numel(a_e2)
        x1=xposa(i);
        y1=yposa(i);
        [ax2,ay2] = find(MeaMap==a_e2(i));
        loc = find(sources==(a1(i)),1,'First');
        score=colors(rank(loc),:);
        sources(loc)=nan; %in case two connections have the same source (non-unique elements in a1)
        al = patchline([x1,xpos(ax2,ay2)],[y1,ypos(ax2,ay2)],'linestyle','-','edgecolor',score,'linewidth',3,'edgealpha',0.7);
%         al = line([x1,xpos(ax2,ay2)],[y1,ypos(ax2,ay2)],'color',score,'linewidth',3);
        uistack(al,'top');
    end
    set(gca,'TickDir','out');
    %--Colorbar Neuros---%
    hn  = axes('visible','off');
    colormap(hn,colorsN);
    c = colorbar('westoutside');
    set(gca,'TickDir','out');
    caxis([floor(min(sN)*10)/10,ceil(max(sN)*10)/10]);
%     ctick = 0:(1/numel(allvalsN))/2:1;
%     ctick=ctick(2:2:end-1);
%     if numel(ctick)>100
% %         set(c,'Ticks',[ctick(1),ctick(floor(end/2)),ctick(end)],'TickLabels',round([allvalsN(1),allvalsN(floor(end/2)),allvalsN(end)].*100)./100);
%         set(c,'Ticks',[0,round(ctick(floor(end/2))*10)/10,max(ceil(ctick))],'TickLabels',[0,round(ctick(floor(end/2))*10)/10,max(ceil(ctick))]);
%     else
%         set(c,'Ticks',ctick(1:2:end),'TickLabels',allvalsN(1:2:end));
%     end
    caxis(caxis);
    set(gca,'FontSize',18);
    
    %--Colorbar Astros---%
    ha  = axes('visible','off');
    colormap(ha,colors);
    c = colorbar;
    set(c,'TickDir','out');
    caxis([floor(min(sA)*10)/10,ceil(max(sA)*10)/10]);
%     ctick = 0:(1/numel(allvals))/2:1;
%     ctick=ctick(2:2:end-1);
%     if numel(ctick)>100
%         set(c,'Ticks',[ctick(1),ctick(floor(end/2)),ctick(end)],'TickLabels',round([allvals(1),allvals(floor(end/2)),allvals(end)].*100)./100);
%     else
%         set(c,'Ticks',ctick(1:2:end),'TickLabels',allvals(1:2:end));
%     end
    caxis(caxis);
    set(gca,'FontSize',18);
end
% set(gcf,'Color','black');
end