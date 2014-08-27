i=5;
traces=DataSet{i}.dfTraces;
time=DataSet{i}.dfTime;
f = savgol(30,5,0);
tr=filtfilt(f,1,traces);
for j=1:size(tr,2)
    DataSet{i}.CaPeaks{j} = peakfinder(tr(:,j),(max(tr(:,j))-min(tr(:,j)))/30);

    for k=1:size(DataSet{i}.CaPeaks{j},1)
        bs(k)  = RollDownBack(tr(:,j),DataSet{i}.CaPeaks{j}(k),mean(tr(:,j))*1.05);
        betry = find(tr(DataSet{i}.CaPeaks{j}(k):end,1)==0,1,'First')+DataSet{i}.CaPeaks{j}(k);
        if isempty(betry)
            be(k)  = RollDownFront(tr(:,j),DataSet{i}.CaPeaks{j}(k),mean(tr(:,j))*1.05);
        else
            be(k)=betry;
        end
    end
%     cutouts=zeros(200,numel(DataSet{i}.PeakOff{j}));
%     for k=1:numel(DataSet{i}.PeakOff{j})
%         cutouts(1:(DataSet{i}.PeakOff{j}(k)-DataSet{i}.PeakOn{j}(k)+1),k)=tr(DataSet{i}.PeakOn{j}(k):DataSet{i}.PeakOff{j}(k),j);
%         cutouts((DataSet{i}.PeakOff{j}(k)-DataSet{i}.PeakOn{j}(k)+2):end,k)=cutouts((DataSet{i}.PeakOff{j}(k)-DataSet{i}.PeakOn{j}(k)+1),k);
%     end
    sames = sort([find(diff(bs)==0),find(diff(bs)==0)+1]);
    [~,loc]=max(tr(sames,1));
    remove = setdiff(sames,sames(loc));
    DataSet{i}.CaPeaks{j}(remove)=[];
    bs(remove)=[];
    be(remove)=[];
    DataSet{i}.PeakOn{j}=bs;
    DataSet{i}.PeakOff{j}=be;
    clear bs; clear be;
end
j=24;
CaRises = time(DataSet{i}.CaPeaks{j});
BurstOnsets = sort([DataSet{i}.bs,DataSet{i}.sbs])./12000;
BurstOffsets = sort([DataSet{i}.be,DataSet{i}.sbe])./12000;
hold on;
plot(time,traces(:,j),'.-');
meanline = ones(size(time)).*mean(traces(:,j));
plot(time,meanline,'--');
plot(time,tr(:,j),'.-');
plot(time(DataSet{i}.CaPeaks{j}),tr(DataSet{i}.CaPeaks{j},j),'.k','markersize',10);
plot(time(DataSet{i}.PeakOn{j}),tr(DataSet{i}.PeakOn{j},j),'.g','markersize',10);
plot(time(DataSet{i}.PeakOff{j}),tr(DataSet{i}.PeakOff{j},j),'.r','markersize',10);

% hold on; plot(time,ones(size(time)).*1.2*mean(tr(:,1)))
% idx = find(traces(:,1)==meanline)
% if isempty(idx)
% idx = find(diff(sign(traces(:,1)-meanline)),1);
% end