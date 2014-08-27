function PlotResultsOnMEA(e1,e2,value,varargin)
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
f = figure;

%---Load Mea Image Data and Map---%
load('MeaMapPlot.mat','MeaImage','MeaMap');
MeaMap = rot90(MeaMap,-1);

%---Plot Mea Image---%
x = repmat([1:size(MeaImage,1)]',1,size(MeaImage,1));
y = repmat(1:size(MeaImage,1),size(MeaImage,2),1);
z = zeros(size(y));
cdata = real2rgb(MeaImage,'gray',[0,255]);
surface(x,y,z,cdata,'EdgeColor','none','FaceColor','texturemap','CDataMapping','direct'); axis('tight');axis('square');

% alpha(0.5);

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
%--Plot Mask--%
if ~isempty(varargin)
    overlay = varargin{1};
    if size(varargin,2)>1
        pos = varargin{2};
        [sitex,sitey]=find(MeaMap==pos);
    else
        sitex=1;
        sitey=1;
    end %must adjust size.
    overlay = imresize(overlay,[357.6,501.2]);
    cdatao = real2rgb(overlay,[255,0,0;152,152,152]./255,[0,1]);
    xo = repmat(1:size(overlay,2),size(overlay,1),1)+xpos(sitex,sitey)-2*size(overlay,1)/3;
    yo = repmat([1:size(overlay,1)]',1,size(overlay,2))+ypos(sitex,sitey)-size(overlay,2)/3-12.5;
    zo = zeros(size(overlay));
    s = surface(xo,yo,zo,cdatao,'EdgeColor','none','FaceColor','texturemap','CDataMapping','direct'); axis('tight');%axis('square');
    xvec = 0:480.3997:2.8824e+03;
    set(gca,'XTick',xvec,'XTickLabel',xvec.*1.0408); % 200/sqrt((xpos(1,3)-xpos(1,2))^2+(ypos(1,3)-ypos(1,2))^2)
    uistack(s,'top');
    alpha(s,'texture');
end

%--Find ROIs and Plot---%
if size(varargin,2)>2
    x1=[]; y1=[]; x2=[]; y2=[];ax2=[];ay2=[];
    a1=varargin{3};
    a_e2=varargin{4};
    avals=varargin{5};
    centers  = cell2mat(arrayfun(@(x) x.Centroid, regionprops(~overlay,'Centroid'),'uniformoutput',0));
    xposa = (centers(a1,1)+xo(1,1))';
    yposa = (centers(a1,2)+yo(1,1))';
    [allvals,rank]=sort([value,avals]);
    [~,~,allv_index] = unique(allvals);
    rank = rank(allv_index);
    if ~isempty(e1) %In case you want to plot just A2N
        a1=a1+max(e1);
        sources = [e1,a1];
    else
        sources = a1;
    end
    sources=sources(rank);
    colors = flip(jet(numel(sources)),1);    
    for i=1:numel(e1)
        [x1,y1] = find(MeaMap==e1(i));
        [x2,y2] = find(MeaMap==e2(i));
        loc = find(sources==e1(i),1,'First');
        score=colors(rank(loc),:);
        sources(loc)=nan;  %in case two connections have the same source (non-unique elements in e1)
        l = line([xpos(x1,y1),xpos(x2,y2)],[ypos(x1,y1),ypos(x2,y2)],'color',score,'linewidth',0.5);
        uistack(l,'top');
    end
    
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
    
    ha  = axes('visible','off');
    colormap(ha,jet(numel(allvals)));
    c = colorbar;
    ctick = 0:(1/numel(allvals))/2:1;
    ctick=ctick(2:2:end-1);
    set(c,'Ticks',ctick,'TickLabels',allvals);
    caxis(caxis);
else
    %---Find Just Electrodes and Plot---%
    colors = jet(numel(e1));
    [value,rank] = sort(value);
    value = round(value*1000)/1000;
    e1=e1(rank);
    e2=e2(rank);
    for i=1:numel(e1)
        [x1,y1] = find(MeaMap==e1(i));
        [x2,y2] = find(MeaMap==e2(i));
        score=colors(rank(i),:);
        l = line([xpos(x1,y1),xpos(x2,y2)],[ypos(x1,y1),ypos(x2,y2)],'color',score,'linewidth',3);
        uistack(l,'top');
    end
    ha  = axes('visible','off');
    colormap(ha,jet(numel(e1)));
    c = colorbar;
    ctick = 0:(1/numel(e1))/2:1;
    ctick=ctick(2:2:end-1);
    set(c,'Ticks',ctick,'TickLabels',value);
    caxis(caxis);
end

end