function DataSet = RestrictToTimeWindow(DataSet,start,stop)
% start  = time in seconds from which to start window
% stop   = time in seconds from which to stop window


winInSamples = [start*12000,stop*12000];
for k=1:size(DataSet,1)
    %-------Cut Time and Time Series-----%
    index = find(DataSet{k}.dfTime>start,1,'First');
    if index>1
        index=index-1;
    end
    
    if ~isempty( find(DataSet{k}.dfTime>stop,1,'First'))
        index(2) = find(DataSet{k}.dfTime>stop,1,'First')-1;
    else
        index(2)=numel(DataSet{k}.dfTime);
    end
    
    DataSet{k}.dfTraces=DataSet{k}.dfTraces(index(1):index(2),:);
    DataSet{k}.FR = DataSet{k}.FR(index(1):index(2),:);
    
    DataSet{k}.dfTime=DataSet{k}.dfTime(index(1):index(2));%-DataSet{k}.dfTime(index(1));
    DataSet{k}.GFR = DataSet{k}.GFR(index(1):index(2));
    
    
    %--------Cut Bursts------------------%
    index(1) = find(DataSet{k}.bs>winInSamples(1),1,'First');
    if index>1
        index=index-1;
    end
    if ~isempty(find(DataSet{k}.bs>winInSamples(2),1,'First'))
        index(2) = find(DataSet{k}.bs>winInSamples(2),1,'First')-1;
    else
        index(2)=numel(DataSet{k}.bs);
    end
    
    
    DataSet{k}.bs=DataSet{k}.bs(index(1):index(2));
    DataSet{k}.be=DataSet{k}.be(index(1):index(2));
    DataSet{k}.bw=DataSet{k}.bw(index(1):index(2));
    
    index(1) = find(DataSet{k}.Trim.bs>winInSamples(1),1,'First');
    if index>1
        index=index-1;
    end
    
    if ~isempty(find(DataSet{k}.Trim.bs>winInSamples(2),1,'First'))
        index(2) = find(DataSet{k}.Trim.bs>winInSamples(2),1,'First')-1;
    else
        index(2)=numel(DataSet{k}.Trim.bs);
    end
    
    
    DataSet{k}.Trim.bs=DataSet{k}.Trim.bs(index(1):index(2));
    DataSet{k}.Trim.be=DataSet{k}.Trim.be(index(1):index(2));
    DataSet{k}.Trim.bw=DataSet{k}.Trim.bw(index(1):index(2));
    
    %---------Cut Spikes------------------%
    [t,ic]=CutSortChannel2(DataSet{k}.t,DataSet{k}.ic,winInSamples(1),winInSamples(2));
    DataSet{k}.t=[];
    DataSet{k}.ic=[];
    DataSet{k}.t=t;
    DataSet{k}.ic=ic;
    [t,ic]=CutSortChannel2(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,winInSamples(1),winInSamples(2));
    DataSet{k}.Trim.t=[];
    DataSet{k}.Trim.ic=[];
    DataSet{k}.Trim.t=t;
    DataSet{k}.Trim.ic=ic;
    
end
end