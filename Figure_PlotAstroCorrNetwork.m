fclose('all');clear all; close all;clc;
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim_noBSinSBS_fixedBursts.mat')
load('CorrDistributions2.mat', 'LagAtMax','MaxCosSim');
cutoff = 0;
dists=[];ins=[]; outs=[];
for which = 1:size(DataSet,1);
    % cutoff = triu(a2a,1)+tril(nan(size(a2a)));
    % cutoff = nanmean(cutoff(:))+nanstd(cutoff(:));
    ic=DataSet{which}.ic;
    a2a=MaxCosSim{which}(size(ic,2)+1:end,size(ic,2)+1:end);
    lagmat = LagAtMax{which}(size(ic,2)+1:end,size(ic,2)+1:end);
    [~, indegree{which},outdegree{which},dist{which},~,~,~,~]=CalcAstroGraph(a2a,lagmat,cutoff);
    temp = dist{which}(~isinf(dist{which}));
    dists = [dists;temp(temp~=0)];
    ins = [ins,indegree{which}];
    outs = [outs,outdegree{which}];
end

%% Plot sample graph
close all;
which=6;
ic=DataSet{which}.ic;
a2a=MaxCosSim{which}(size(ic,2)+1:end,size(ic,2)+1:end);
lagmat = LagAtMax{which}(size(ic,2)+1:end,size(ic,2)+1:end);
cutoff = triu(a2a,1)+tril(nan(size(a2a)));
cutoff = nanmean(cutoff(:))+1.5*nanstd(cutoff(:));
[bg, ~,~,~,~,emptynodes,noderank,edgerank]=CalcAstroGraph(a2a,lagmat,0);
view(bg);
temp = a2a;
temp(temp<cutoff)=0;
PlotBGrealspace(DataSet{which}.mask,temp,emptynodes,noderank);
%% Plot shortest path, in and out degree
figure;
subplot = @(m,n,p) subtightplot (m, n, p, [0.07 0.05], [0.05 0.05], [0.05 0.05]);
subplot(1,2,1)
% dbins = OPTBINS(dists',500);
dbins = 458;
dispHIST(dists',dbins);
xlabel('weighted directed shortest path');
axis tight;
xlim([0,1]);
set(gca,'TickDir','out','PlotBoxAspectRatio',[1,1,1]);
box off;
subplot(1,2,2)
% ibins = OPTBINS(ins,500);
ibins = 6;
hold on;
[h1,e,c] = PlotHistWithErrors(ins,ibins,1,'b',[],[]);
[h2,~] = PlotHistWithErrors(outs,ibins,0.75,'r',e,c);
% y  = pdf(PD{1}.DistributionName,linspace(bins(1),bins(end),numel(bins)*10),PD{1}.mu,PD{1}.omega);
% h3 = plot(linspace(bins(1),bins(end),numel(bins)*10),y,'k','LineWidth',3);
axis tight;
xlabel('weighted degree');
ylabel('probability density');
% legend([h1,h2,h3],{'in degree','out degree','nagakami fit'});
legend([h1,h2],{'in degree','out degree'});
set(gca,'TickDir','out','FontSize',9,'XTick',c,'XTickLabel',round(c*100)/100,'PlotBoxAspectRatio',[1,1,1]);
box off;
%% Plot 2D Histogram of Correlation vs. Distance
figure;
load('CorrDist2.mat','a2as','a2ad');
s = cell2mat(a2as');
d = cell2mat(a2ad');
optM = [15,13];
bin = {linspace(min(d),max(d),optM(1)),linspace(min(s),max(s),optM(2))};
normhist=zeros(numel(bin{1}),numel(bin{2}));
for i=1:size(a2as,2)
    [temphist,bins]=hist3([a2ad{i},a2as{i}],'Edges',bin);
    distmat = regiondist(~DataSet{i}.mask);
    
    
    ix = repmat(1:size(distmat,1),[size(distmat,1),1]); %Create index
    ix = triu(ix);
    ix = ix+tril(nan(size(ix)));
    ix =ix(~isnan(ix));
    
    factor=zeros(size(distmat,1),optM(1));
    bincheck = [bins{1},inf];
    for j=1:size(distmat,1)
        for k=1:optM(1)
            factor(j,k) =  sum((distmat(j,:)>=bincheck(k))&(distmat(j,:)<bincheck(k+1)));
        end
    end
    
    fac = zeros(optM(1),optM(2));
    bincheck2 = [bins{2},inf];
    for x = 1:optM(1)
        set1 = ix(find(a2ad{i}>=bincheck(x) & a2ad{i}<bincheck(x+1))); %distance bin
        for y=1:optM(2)
            set2 = ix(find(a2as{i}>=bincheck2(y) & a2as{i}<bincheck2(y+1))); %correlation score bin
            [~,~,include] = intersect(sort(set1),sort(set2));
            include = set2(include);
            fac(x,y)=sum(factor(include,x),1);
        end
    end
    
    
    fac(find(fac==0))=inf;
    temphist = temphist./fac;
    normhist=normhist+temphist/trapz(bins{2},trapz(bins{1},temphist));
end
normhist = normhist./size(a2as,2);
xo = repmat(bins{1}',1,size(normhist,2));
yo = repmat(bins{2},size(normhist,1),1);
surface(xo,yo,normhist);
set(gca,'PlotBoxAspectRatio',[optM(1),optM(2),1],'TickDir','out');
az=0;%35
el=90;%32
view([az,el]);
% set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
ylabel('normalized unbiased astrocyte cross-correlation');
xlabel('distance [um]');
grid off;
axis tight;
set(gca,'FontSize',9,'TickDir','out');
c= colorbar;
ylabel(c,'probability density');
set(c,'FontSize',9,'TickDir','Out');
