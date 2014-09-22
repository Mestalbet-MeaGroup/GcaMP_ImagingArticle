screen_size = get(0, 'ScreenSize');
set(gcf, 'PaperUnits', 'centimeters');
for i=3%:9
    tr = DataSet{i}.dfTraces;
        for j=1:size(tr,2)
          p = ConservativePeakDetection(tr(:,j));
          figure('visible','off');
          if ~isempty(p)
            hold on; plot(tr(:,j)); plot(p,tr(p,j),'or');
            set(gcf, 'PaperPosition', [0 0 15 10]); %x_width=10cm y_width=15cm
            saveas(gcf,['TraceDetect_Cult_',num2str(DataSet{i}.culture),' ','Trace_',num2str(j),'.tif']);
%             export_fig(['TraceDetect_Cult_',num2str(DataSet{i}.culture),' ','Trace_',num2str(j),'.tif'],'-r300');
          end
        close all;  
        end
end
