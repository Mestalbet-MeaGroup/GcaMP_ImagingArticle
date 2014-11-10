function DataSet = FixBursts(DataSet);

for j=5:size(DataSet,1)
    t= DataSet{j}.t;
    ic = DataSet{j}.ic;
    bs = DataSet{j}.bs;
    sbs = DataSet{j}.sbs;
    sbe = DataSet{j}.sbe;
  
    
    if ~isempty(sbs)
        [tt,tic] = RemoveHAneurons(t,ic,[bs,sbs],[bs+4000,sbe]);
        [bs,be,~,~,~,~]=UnsupervisedBurstDetection2(tt,tic);
        for i=1:numel(sbs)
            places = find(be>sbs(i) & bs<=sbe(i));
            bs(places)=[];
            be(places)=[];
        end
    else
        [tt,tic] = RemoveHAneurons(t,ic,bs,bs+4000);
        [bs,be,~,~,~,~]=UnsupervisedBurstDetection2(tt,tic);
    end
    newbs=[];
    newbe=[];
    vec=ConvertIC2Samora(tic);
    for i=1:numel(bs)
        [temp,ix] = sort(tt);
        vec=vec(ix);
        locs = find(temp>=bs(i) & temp<=be(i));
        numchannels = numel(unique(vec(locs)));
        if numchannels>=(0.5*size(ic,2))
            newbs = [newbs,bs(i)];
            newbe = [newbe,be(i)];
        end
    end
    DataSet{j}.bs=newbs;
    DataSet{j}.be=newbe;
    DataSet{j}.bw=newbe-newbs;
end
end