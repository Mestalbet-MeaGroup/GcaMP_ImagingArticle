for i=1:size(DataSet,1)
    PlotAstroTracesWithRaster(DataSet{i}.dfTraces',DataSet{i}.dfTime,DataSet{i}.t,DataSet{i}.ic);
    export_fig(['RasterTrace_' num2str(i) '.tif'],'-native');
    close all;
end

% for i=1:size(DataSet,1)
%     J = entropyfilt(DataSet{i}.dfTraces);
%     [a,b]=CalcDf_f(J',DataSet{i}.fs,DataSet{i}.dfTime);
%     PlotAstroTracesWithRaster(a,b,DataSet{i}.t,DataSet{i}.ic);
%     export_fig(['RasterTraceEntropy_' num2str(i) '.tif'],'-native');
%     close all;
% end

% traces = DataSet{i}.RawTraces;
% traces(:,26303)=traces(:,26302);
% [dfTraces,dfTime]=CalcDf_f(traces,DataSet{i}.fs,DataSet{i}.RawTime);
J = entropyfilt(dfTraces');
[a,b]=CalcDf_f(J',DataSet{i}.fs,dfTime);
PlotAstroTracesWithRaster(dfTraces',dfTime,DataSet{i}.t,DataSet{i}.ic);
% fun = @(block_struct) block_struct.data-mean(block_struct.data);
% I2 = blockproc(dfTraces',[100 1],fun);