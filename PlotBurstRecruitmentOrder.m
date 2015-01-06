function PlotBurstRecruitmentOrder(channelOrderR,channelOrderNR)
%CREATEFIT    Create plot of datasets and fits
%   CREATEFIT(CHANNELORDERR,CHANNELORDERNR)
%   Creates a plot, similar to the plot in the main distribution fitting
%   window, using the data that you provide as input.  You can
%   apply this function to the same data you used with dfittool
%   or with different data.  You may want to edit the function to
%   customize the code and this help message.
%
%   Number of datasets:  2
%   Number of fits:  0
%
%   See also FITDIST.

% This function was automatically generated on 11-Dec-2014 13:53:53

% Data from dataset "channelOrderR data":
%    Y = channelOrderR

% Data from dataset "channelOrderNR data":
%    Y = channelOrderNR

% Force all inputs to be column vectors
channelOrderR = channelOrderR(:);
channelOrderNR = channelOrderNR(:);

% Prepare figure
clf;
hold on;
LegHandles = []; LegText = {};


% --- Plot data originally in dataset "channelOrderR data"
[CdfF,CdfX] = ecdf(channelOrderR,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 1;
[~,BinEdge] = internal.stats.histbins(channelOrderR,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
    'LineStyle','-', 'LineWidth',1);
xlabel('Data');
ylabel('Density')
LegHandles(end+1) = hLine;
LegText{end+1} = 'Spike order: bursts with Ca';

% --- Plot data originally in dataset "channelOrderNR data"
[CdfF,CdfX] = ecdf(channelOrderNR,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 1;
[~,BinEdge] = internal.stats.histbins(channelOrderNR,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'FaceColor','none','EdgeColor',[0.333333 0.666667 0],...
    'LineStyle','-', 'LineWidth',1);
xlabel('Rank orders of firing for view-field electrodes in a burst');
ylabel('Density')
LegHandles(end+1) = hLine;
LegText{end+1} = 'Spike order: bursts without Ca'

% Create grid where function will be computed
XLim = get(gca,'XLim');
XLim = XLim + [-1 1] * 0.01 * diff(XLim);
XGrid = linspace(XLim(1),XLim(2),100);


% Adjust figure
box on;
hold off;

% Create legend from accumulated handles and labels
hLegend = legend(LegHandles,LegText,'Orientation', 'vertical', 'FontSize', 9, 'Location', 'northeast');
set(hLegend,'Interpreter','none');
