function DataSet = CalcTrimTIC(DataSet)
h = waitbarwithtime(0,'Please Wait...');
for i=1:size(DataSet,1)
    [t,ic,DataSet{i}.Trim.bs,DataSet{i}.Trim.be,DataSet{i}.Trim.bw] = SortChannelsByFR2(DataSet{i}.t,DataSet{i}.ic,DataSet{i}.bs,DataSet{i}.be,DataSet{i}.bw);
    for j=1:numel(DataSet{i}.sbs)
        delete = find( (DataSet{i}.Trim.bs>=DataSet{i}.sbs(j)) & (DataSet{i}.Trim.bs<DataSet{i}.sbe(j)) );
        DataSet{i}.Trim.bs(delete)=[];
        DataSet{i}.Trim.be(delete)=[];
        DataSet{i}.Trim.bw(delete)=[];
    end
    [DataSet{i}.Trim.t,DataSet{i}.Trim.ic] = OnlyNonHA(t,ic,[DataSet{i}.Trim.bs,DataSet{i}.sbs],[DataSet{i}.Trim.be,DataSet{i}.sbe]);
    waitbarwithtime(i/size(DataSet,1),h);
end

%% Burst Detection Verification
% i=1;
% j=4;
% % PlotRasterWithBursts(DataSet{j}.Trim.t,DataSet{j}.Trim.ic,DataSet{j}.Trim.bs,DataSet{j}.Trim.be);
% PlotRasterWithBursts(DataSetStims{j}.Trim.t,DataSetStims{j}.Trim.ic,DataSetStims{j}.Trim.bs,DataSetStims{j}.Trim.be);
end