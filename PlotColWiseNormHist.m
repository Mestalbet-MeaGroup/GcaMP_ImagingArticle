function [normhist,bins] = PlotColWiseNormHist(data1,data2,varargin);
%inputs: two data vectors with which to calculate the 2D histogram
%option inputs: 1/ 2D vector number of bins [bin1,bin2]
%               2/ normalize counts by col, row or both
subplot = @(m,n,p) subtightplot (m, n, p, [0.005 0.005], [0.05 0.05], [0.05 0.05]);

if isempty(varargin)
    optM = nsOPTBINS([cell2mat(fr);cell2mat(Pca)]');
else
    optM=varargin{1};
end

if size(data1,1)>size(data1,2)
    data1=data1';
    if size(data1,1)~=size(data2,1)
        data2=data2';
    end
end

[temphist,bins]=hist3([data1;data2]',optM([2,1]));
hc = sum(temphist,2);
hc2 = sum(temphist,1);

temphist1=temphist;

if numel(varargin)>1
    switch varargin{2}
        case 'both'
            %---Normalize by total row counts---%
%             counts = repmat(sum(temphist,2),[1,optM(1)]);
%             counts(counts==0)=inf;
%             temphist1 = temphist./counts;
%             %----------------------------%
%             
%             %---Normalize by total column counts---%
%             counts2 = repmat(sum(temphist,1),[optM(2),1]);
%             counts2(counts2==0)=inf;
%             temphist = temphist1./counts2;
            %----------------------------%
            
            counts = repmat(sum(temphist,2),[1,optM(1)]);
            counts = counts/sum(counts(:,1));
            counts2 = repmat(sum(temphist,1),[optM(2),1]);
            counts2 = counts2/sum(counts2(1,:));
            expected = counts.*counts2;
            expected(expected==0)=inf;
            temphist = temphist./expected;
            
        case 'col'
            %---Normalize by total row counts---%
            counts = repmat(sum(temphist,2),[1,optM(1)]);
            counts = counts / sum(counts(:,1));
            counts(counts==0)=inf;
            temphist1 = temphist./counts;
            %----------------------------%
        case 'row'
            %---Normalize by total column counts---%
            counts2 = repmat(sum(temphist,1),[optM(2),1]);
            counts2 = counts2/sum(counts2(1,:));
            counts2(counts2==0)=inf;
            temphist = temphist1./counts2;
            %----------------------------%
        case 'none'
            display('no count normalization');
        otherwise
            error('Invalid Normalization Selection');
    end
end

subplot(6,6,2:6)
% [hc,bc]=histcounts(data1,bins{1});
% if mod(numel(bc),2)==0
%     bar(bc,[hc,0],'histc');
% else
%     bar(bc,[0,hc],'histc');
% end
stairs(bins{1},hc,'LineWidth',2);
axis tight;
xlimits = get(gca,'XLim');
% set(allchild(gca),'FaceColor',[0 0 0]);
set(gca,'XTick',[],'YTick',[]);

subplot(6,6,[7,13,19,25,31])
stairs(bins{2},hc2,'LineWidth',2);
axis tight;
ylimits = get(gca,'XLim');
set(gca,'YTick',[],'YDir','reverse','XDir','reverse','TickDir','out','TickLength',[0.005,0.005]);
view([90,90]);
% set(allchild(gca),'FaceColor',[0 0 0]);
set(gca,'FontSize',9,'ycolor','w');

subplot(6,6,[8:12,14:18,20:24,26:30,32:36]);
%---Create Surface---%
normhist = temphist/trapz(bins{2},trapz(bins{1},temphist));
xo = repmat(bins{1}',1,size(normhist,2));
yo = repmat(bins{2},size(normhist,1),1);
if numel(varargin)>2
    if strcmp(varargin{3},'log')
        if min(normhist(normhist>0))<1
            factor = ceil(abs(log10(min(normhist(normhist>0))))); %for when values are less than one
        else
            factor=0;
        end
        surface(xo,yo,normhist*10^factor);
        xlim(xlimits);
        ylim(ylimits);
        az=0;%35
        el=90;%32
        view([az,el]);
        set(get(gca,'child'),'FaceColor','flat','CDataMode','auto');
        set(gca,'YTick',[]);
        CData = get(get(gca,'child'),'CData');
        minC = nanmin(CData(CData>0)); % Find min
        maxC = nanmax(CData(:)); % Find max
        CData = log10(CData); % Convert data to log space
        Data(isinf(CData))=0;
        % Now set the CData of the surface to the normalized CData
        set(get(gca,'child'),'CData',CData,'CDataMapping','scaled')
        grid off;
        set(gca,'TickDir','out','FontSize',9);
        hCbar = colorbar(gca,'eastoutside');
        range = [-inf,0:factor];
        ylim(hCbar,[0,factor]);
        if min(normhist(normhist>0))<1
            ticks = get(hCbar,'YTick').*-1;
            ticks(1)=0;
            ticks = ticks([end:-1:1]);
            set(hCbar,'YTickLabel',ticks,'TickDir','out');
        else
            set(hCbar,'YTickLabel',get(gca,'YTick'),'TickDir','out');
        end
    end
else
    surface(xo,yo,normhist);
    xlim(xlimits);
    ylim(ylimits);
    az=0;%35
    el=90;%32
    view([az,el]);
    set(get(gca,'child'),'FaceColor','flat','CDataMode','auto');
    set(gca,'YTick',[],'TickDir','out','FontSize',9);
    colormap([[0,0,0];parula(numel(unique(normhist(:))))]);
    hCbar = colorbar(gca,'eastoutside');
end


