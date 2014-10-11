function [e1,e2,eval,a1,a_e2,aval,ch,mask]=ParseConnectionValuesWithVF(which,number,varargin);
%Function which allows you to choose which culture and how many top correlation
%values. Used to produce input variables for PlotResultsOnMea.m. If you
%want to use another value than maximum correlation (from
%CorrDistributions), add the matrix as an optional third input.
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
load('CorrDistributions2.mat', 'MaxCosSim');
% load('CorrDistribution_Cult6_temp2.mat','MaxCosSim');
MaxCosSim = MaxCosSim{which};
numElec = size(DataSet{which}.ic,2);
load('MeaMapPlot.mat','MeaMap');
ic =DataSet{which}.ic(1,:);
%--N2N---%
if isempty(varargin)
    temp = MaxCosSim(1:numElec,1:numElec);
    temp = triu(temp,1);
else
    temp = triu(varargin{1},1);
end

for i=1:number
    [eval(i),loc(i)]=nanmax(temp(:));
    temp(loc(i))=nan;
end

if isempty(varargin)
    temp = MaxCosSim(1:numElec,1:numElec);
    temp = triu(temp,1);
    temp = temp+temp';
else
    temp = triu(varargin{1},1);
    temp = temp+temp';
end

[~,~,channels] = FindElecsinVF(DataSet{which}.channel,MeaMap);
[a,b]=ind2sub(size(temp),loc);
e1 = DataSet{which}.ic(1,a);
e2 = DataSet{which}.ic(1,b);
channels(ismember(channels,e1))=[];
channels(ismember(channels,e2))=[];
if ~isempty(channels)
    vfelec = temp(ismember(ic,channels),:);
    channels = channels(ismember(channels,ic));
    [c,d]=find(vfelec>0);
    [~,c]=find(bsxfun(@eq, channels(c),ic));
    e1= [e1,ic(c)];
    e2= [e2,ic(d)];
    temp2 = arrayfun(@(x,y) temp(x,y),c,d); 
    eval= [eval,temp2'];
end

%--A2N--%
temp = MaxCosSim(1:numElec,numElec+1:end);
for i=1:number
    [aval(i),loc(i)]=nanmax(temp(:));
    temp(loc(i))=nan;
end
[a,a1]=ind2sub(size(temp),loc);
a_e2=DataSet{which}.ic(1,a);
% [~,~,channels] = FindElecsinVF(DataSet{which}.channel,MeaMap);
% channels(ismember(channels,a_e2))=[];
% if ~isempty(channels)
%     vfelec = temp(ismember(DataSet{which}.ic(1,:),channels),:);
%     [c,d]=find(vfelec>0);
%     [~,c]=find(bsxfun(@eq, channels(c),DataSet{which}.ic(1,:)));
%     a_e2= [a_e2,DataSet{which}.ic(1,c)];
%     a1= [a1,d'];
%     temp2 = arrayfun(@(x,y) temp(x,y),c,d); 
%     aval= [aval,temp2'];
% end
ch = DataSet{which}.channel;
mask = DataSet{which}.mask;
end


