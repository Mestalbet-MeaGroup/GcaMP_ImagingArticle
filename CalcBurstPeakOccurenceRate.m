function [numBI,numNoB,numSB]=CalcBurstPeakOccurenceRate(DataSet,k);

[b2pi,p2bi,p2b_cv,b2p_cv,peaks]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,DataSet{k}.bs,sort(DataSet{k}.t));
if ~isempty(DataSet{k}.sbs)
    [sb2pi,p2sbi,p2sb_cv,sb2p_cv,peaksSB]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,DataSet{k}.sbs,sort(DataSet{k}.t));
end
%% Find Set of Peaks with Bursts Nearby
for i=1:size(p2bi,2)
    dis = p2b_cv{i};
    index = find((dis>0)&(dis<1));
    bi{i} = p2bi{i}(index);
    pi{i} = peaks{i}(index);
end

%% Find Set of Peaks with Bursts Far Away
if ~isempty(DataSet{k}.sbs)
    for i=1:size(p2bi,2)
        dis = p2b_cv{i};
        index1 = find(abs(dis)>3);
        dis2 = abs(p2sb_cv{i});
        sbw = DataSet{k}.sbw(p2sbi{i})./12000;
        index2 = find(dis2<sbw);
        index =index1(~ismember(index1,index2));
        Nob_pi{i} = peaks{i}(index);
    end
else
    for i=1:size(p2bi,2)
        dis = p2b_cv{i};
        index1 = find(abs(dis)>3);
        Nob_pi{i} = peaks{i}(index1);
    end
end

%% Find Set of Bursts with no corresponding Peaks
if ~isempty(DataSet{k}.sbs)
    for i=1:size(p2bi,2)
        dis = b2p_cv{i};
        index1 = find(abs(dis)>3);
        dis2 = abs(sb2p_cv{i});
        sbw = DataSet{k}.sbw(p2sbi{i}(sb2pi{i}))./12000;
        index2 = find(dis2<sbw);
        index =index1(~ismember(index1,index2));
        No_bi{i} = index;
    end
else
    for i=1:size(p2bi,2)
        dis = b2p_cv{i};
        index1 = find(abs(dis)>3);
        %     dis2 = abs(sb2p_cv{i});
        %     sbw = DataSet{k}.sbw(p2sbi{i}(sb2pi{i}))./12000;
        %     index2 = find(dis2<sbw);
        %     index =index1(~ismember(index1,index2));
        No_bi{i} = index1;
    end
end
%% Find Set of Peaks near SuperBursts
if ~isempty(DataSet{k}.sbs)
    
    for i=1:size(p2sbi,2)
        dis = p2sb_cv{i};
        index = find((dis>0)&(dis<3));
        SB_pi{i} =  peaks{i}(index);
        SB_bi{i} =  unique(p2sbi{i}(index));
    end
    
    numBI = numel(unique(cell2mat(bi)))./numel(DataSet{k}.bs); % Percentage of bursts with a corresponding calcium increase
    temp=[];
    for i=1:size(Nob_pi,2),
        if ~isempty(Nob_pi{i}),
            temp = [temp;unique(Nob_pi{i})];
        end
    end
end
%%
numBI = numel(unique(cell2mat(bi)))./numel(DataSet{k}.bs); % Percentage of bursts with a corresponding calcium increase
temp=[];
for i=1:size(Nob_pi,2),
    if ~isempty(Nob_pi{i}),
        temp = [temp;unique(Nob_pi{i})];
    end
end
numNoB = numel(temp)./numel(cell2mat(cellfun(@(x) unique(x),peaks,'UniformOutput',0)')); % Percentage of calcium increases without a corresponding burst
c=1;
if ~isempty(DataSet{k}.sbs)
    for i=1:size(SB_bi,2)
        if ~isempty(SB_bi{i})
            n{c}=SB_bi{i};
            c=c+1;
        end
    end
    try
        numSB = numel(unique(cell2mat(n)))./numel(DataSet{k}.sbs); % Percentage of superbursts with a corresponding calcium increase
    catch
        numSB = numel(unique(cell2mat(n')))./numel(DataSet{k}.sbs); % Percentage of superbursts with a corresponding calcium increase
    end
else
    numSB=nan;
end
end