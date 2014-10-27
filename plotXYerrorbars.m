function [h]=plotXYerrorbars(x, y, ex, ey,varargin)

% create the standard  scatterplot
h=plot(x, y, varargin{:});

% looks better with large points
set(h, 'MarkerSize', 25);
hold on

% now draw the vertical errorbar for each point
for i=1:length(x)
        xV = [x(i); x(i)];
        yV = [y(i); y(i)];
      
        xMin = x(i) + ex(i);
        xMax = x(i) - ex(i);
        yMin = y(i) + ey(i);
        yMax = y(i) - ey(i);
     
        xB = [xMin, xMax];
        yB = [yMin, yMax];
         % draw error bars
%         h=plot(xV, yV, varargin{:});
%         set(h, 'LineWidth', 1);
        h=plot(xB, yV,varargin{:});
        set(h, 'LineWidth', 1);
        h=plot(xV, yB,varargin{:});
        set(h, 'LineWidth', 1);
        
        
end
    ha = patch([x fliplr(x)], [ex,ey],'r');