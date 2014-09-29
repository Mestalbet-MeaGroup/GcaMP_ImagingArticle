% a2a= tril(nan(size(a2a)))+triu(a2a,1);
close all; clear all;
loadcult(6);
dg=eye(size(a2a));
dg(dg==1)=nan;
ameans=nanmean(a2a+a2a'+dg,1); 
amaxs=nanmax(a2a+a2a'+dg,[],1); 
areas = arrayfun(@(x) x.Area,regionprops(~mask,'Area'))';
acents=pagerank(triu(a2a,1)+triu(a2a,1)');

figure;
[~,h1,h2]=plotyy(areas,ameans,areas,acents);
set(h1,'LineStyle','none','Marker','o');
set(h2,'LineStyle','none','Marker','+');
title('Areas vs. Means/Centrality');
legend({'Mean','Centrality'});

figure;
[~,h1,h2]=plotyy(areas,amaxs,areas,acents);
set(h1,'LineStyle','none','Marker','o');
set(h2,'LineStyle','none','Marker','+');
title('Areas vs. Max/Centrality');
legend({'Max','Centrality'});

figure;
[~,h1,h2]=plotyy(acents,ameans,acents,amaxs);
set(h1,'LineStyle','none','Marker','o');
set(h2,'LineStyle','none','Marker','+');
title('Centrality vs. Mean/Max');
legend({'Mean','Max'});
