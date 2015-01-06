function PlotBurstRecruitment(ratioR,ratioNR)
%CREATEFIT    Create plot of datasets and fits
%   CREATEFIT(RATIOR,RATIONR)
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

% This function was automatically generated on 12-Dec-2014 12:18:26

% Data from dataset "ratioR data":
%    Y = ratioR

% Data from dataset "ratioNR data":
%    Y = ratioNR

% Force all inputs to be column vectors
ratioR = ratioR(:);
ratioNR = ratioNR(:);

% Prepare figure
clf;
hold on;
LegHandles = []; LegText = {};


% --- Plot data originally in dataset "ratioR data"
[CdfF,CdfX] = ecdf(ratioR,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 1;
[~,BinEdge] = internal.stats.histbins(ratioR,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'FaceColor','none','EdgeColor',[0.333333 0 0.666667],...
    'LineStyle','-', 'LineWidth',1);
xlabel('Data');
ylabel('Density')
LegHandles(end+1) = hLine;
LegText{end+1} = 'Bursts with Ca oscillations';

% --- Plot data originally in dataset "ratioNR data"
[CdfF,CdfX] = ecdf(ratioNR,'Function','cdf');  % compute empirical cdf
BinInfo.rule = 1;
[~,BinEdge] = internal.stats.histbins(ratioNR,[],[],BinInfo,CdfF,CdfX);
[BinHeight,BinCenter] = ecdfhist(CdfF,CdfX,'edges',BinEdge);
hLine = bar(BinCenter,BinHeight,'hist');
set(hLine,'FaceColor','none','EdgeColor',[0.333333 0.666667 0],...
    'LineStyle','-', 'LineWidth',1);
xlabel('Percentage view-field participation in bursts');
ylabel('Density')
LegHandles(end+1) = hLine;
LegText{end+1} = 'Bursts without Ca oscillations';

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
