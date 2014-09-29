figure;
hold on;
vec = [nanmean(n2n,2),nanmean(a2n,2)];
vec(isnan(vec))=0;
hist3(vec,nsOPTBINS(vec'));
n = hist3(vec,nsOPTBINS(vec'));
n1 = n';
n1(size(n,2) + 1, size(n,1) + 1) = 0;
xb = linspace(min(vec(:,1)),max(vec(:,1)),size(n,1)+1);
yb = linspace(min(vec(:,2)),max(vec(:,2)),size(n,2)+1);
h = pcolor(xb,yb,n1);
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot)
grid on
view(3);
axis tight;
set(gcf,'Color','w');
set(gca,'FontSize',18);
xlabel('Neuron Correlation in Neuron Network');
ylabel('Neuron Correlation in Astro-Neuro Network');

figure;
hold on;
vec = [nanmean(a2a,2),nanmean(a2n,1)'];
vec(isnan(vec))=0;
vec(isnan(vec))=0;
hist3(vec,nsOPTBINS(vec'));
n = hist3(vec,nsOPTBINS(vec'));
n1 = n';
n1(size(n,2) + 1, size(n,1) + 1) = 0;
xb = linspace(min(vec(:,1)),max(vec(:,1)),size(n,1)+1);
yb = linspace(min(vec(:,2)),max(vec(:,2)),size(n,2)+1);
h = pcolor(xb,yb,n1);
h.ZData = ones(size(n1)) * -max(max(n));
colormap(hot)
grid on
view(3);
axis tight;
set(gcf,'Color','w');
set(gca,'FontSize',18);
% xlabel('Astrocyte Correlation in Astro Network');
% ylabel('Astrocyte Correlation in Astro-Neuro Network');