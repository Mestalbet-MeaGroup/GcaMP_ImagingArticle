
fun1  = @(x) triu(x,1)+tril(nan(size(x)),1);
fun2 = @(x) x(:);

corvals = unique(cell2mat(cellfun(fun2,cellfun(fun1, MaxCosSim,'UniformOutput',false),'UniformOut',false)));
corvals = corvals(~isnan(corvals));
corvals = corvals(corvals>0.01);
corvals = corvals(corvals<1);

[~,bins]=hist(corvals,100);
for i=1:9
    corrs = MaxCosSim{i};
    corrs = triu(corrs,1)+tril(nan(size(corrs)),1);
    numc = size(DataSet{i}.ic,2);
    n2n = corrs(1:numc,1:numc);
    n2n = n2n(:);
    n2n = n2n(~isnan(n2n));
    n2nh(:,i)=histc(n2n,bins);
    n2nh(:,i)=n2nh(:,i)./sum(n2nh(:,i));
    
    a2a = corrs(numc+1:end,numc+1:end);
    a2a = a2a(:);
    a2a = a2a(~isnan(a2a));
    a2ah(:,i)=histc(a2a,bins);
    a2ah(:,i)=a2ah(:,i)./sum(a2ah(:,i));
    
    a2n = corrs(1:numc,numc+1:end);
    a2n = a2n(:);
    a2n = a2n(~isnan(a2n));
    a2nh(:,i)=histc(a2n,bins);
    a2nh(:,i)=a2nh(:,i)./sum(a2nh(:,i));
end
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.005], [0.05 0.05], [0.05 0.05]);
subplot(3,1,1)
hold on;
h=bar(bins,nanmean(n2nh,2),'histc');
plot(bins,(n2nh./max(n2nh,2))*.5,'--');
set(gca,'XTick',[]);
set(gca,'FontSize',18)
axis('tight');
xlim([0,1]);

subplot(3,1,2)
hold on;
h=bar(bins,nanmean(a2nh,2),'histc');
x = get(h,'XData');
plot(bins,(n2nh./max(a2nh,2))*1.1,'--');
set(gca,'FontSize',18)
set(gca,'XTick',[]);
axis('tight');
xlim([0,1]);

subplot(3,1,3)
hold on;
h=bar(bins,nanmean(a2ah,2),'histc');
plot(bins,(n2nh./max(a2ah,2))*.5,'--');
set(gca,'FontSize',18)
axis('tight');
xlim([0,1]);
%%
nnscores = [];
nascores = [];
anscores = [];
aascores = [];
for i=1:9
    corrs = MaxCosSim{i};
    corrs = triu(corrs,1)+tril(nan(size(corrs)),1);
    numc = size(DataSet{i}.ic,2);
    n2n = corrs(1:numc,1:numc);
    a2a = corrs(numc+1:end,numc+1:end);
    a2n = corrs(1:numc,numc+1:end);
    
%     nnscores = [nnscores;nanmean(n2n,2)];
%     nascores = [nascores;nanmean(a2n,2)];
%     anscores = [anscores,nanmean(a2n,1)];
%     aascores = [aascores;nanmean(a2a,2)];
%     
    %     for j=1:size(n2n,1), numspks(j) = numel(DataSet{i}.t(DataSet{i}.ic(3,j):DataSet{i}.ic(4,j))); end
        nnscores = [nnscores;max(n2n,[],2)];
        nascores = [nascores;max(a2n,[],2)];
        anscores = [anscores,max(a2n,[],1)];
        aascores = [aascores;max(a2a,[],2)];
    
    
end
steps=0.05;
window=0.1;

for i=0:floor(max(nnscores)/steps)

    start  = 0+i*steps;
    stop = start+window;
    nns(i+1) = mean( nnscores( (nnscores>=start) & (nnscores<=stop) ) );
    Enns(i+1) = std( nnscores( (nnscores>=start) & (nnscores<=stop) ) );
    asns(i+1) = mean( nascores( (nnscores>=start) & (nnscores<=stop) ) );
    Easns(i+1) = std( nascores( (nnscores>=start) & (nnscores<=stop) ) );
end


for i=0:floor(max(aascores)/steps)

    start  = 0+i*steps;
    stop = start+window;
    aas(i+1) = mean( aascores( (aascores>=start) & (aascores<=stop) ) );
    Eaas(i+1) = std( aascores( (aascores>=start) & (aascores<=stop) ) );
    nas(i+1) = mean( anscores( (aascores>=start) & (aascores<=stop) ) );
    Enas(i+1) = mean( anscores( (aascores>=start) & (aascores<=stop) ) );
end
% hold on;
% h=plotXYerrorbars(nns, asns, Enns, Easns,'r');
% h=plotXYerrorbars(aas, nas, Eaas, Enas,'k');
% 

% plot(nns,ans,'.-');
% plot(aas,nas,'.-');
[~,bins] = hist3([nnscores,nascores],[50,25]);
subplot(1,2,1);
hist3([nnscores,nascores],'edges',bins);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
axis tight;
view([0,90]);
ylim([min(bins{2}),max(bins{2})]);
xlim([min(bins{1}),max(bins{1})]);
set(gca,'FontSize',18);

subplot(1,2,2);
hist3([aascores,anscores'],'edges',bins);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
axis tight;
view([0,90]);
set(gca,'YTick',[]);
ylim([min(bins{2}),max(bins{2})]);
xlim([min(bins{1}),max(bins{1})]);
set(gca,'FontSize',18);
% xlim([0,1]);
% ylim([0,0.7]);
