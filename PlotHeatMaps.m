load('CorrDist.mat');
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat');
subplot = @(m,n,p) subtightplot (m, n, p, [0.05 0.05], [0.05 0.05], [0.05 0.05]);

%---Correlation vs. Distance---%
figure;
xlimit = 500;
for i=1:9
    subplot(3,3,i)
    optM = a2aB{i};
    a2a = [a2ad{i}(~isnan(a2as{i})),a2as{i}(~isnan(a2as{i}))]';
    a2a = a2a(:,a2a(2,:)~=0);
    Plot2DHist(a2a(1,:)',a2a(2,:)',optM(1),optM(2),xlimit);
    title(['A2A: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end

figure;
xlimit = 4000;
for i=1:9
    subplot(3,3,i)
    optM = a2nB{i};
    a2n = [a2nd{i}(~isnan(a2ns{i})),a2ns{i}(~isnan(a2ns{i}))]';
    a2n = a2n(:,a2n(2,:)~=0);
    Plot2DHist(a2n(1,:)',a2n(2,:)',optM(1),optM(2),xlimit);
    title(['A2N: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end

figure;
for i=1:9
    subplot(3,3,i)
    optM = n2nB{1};
    n2n = [n2nd{i}(~isnan(n2ns{i})),n2ns{i}(~isnan(n2ns{i}))]';
    n2n = n2n(:,n2n(2,:)~=0);
    Plot2DHist(n2n(1,:)',n2n(2,:)',round(optM(1)/3),optM(2)*2,xlimit);
    title(['N2N: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end
clear_all_but('DataSet');
load('LagDist.mat');

%---Lag versus Distance---%
figure;
xlimit = 500;
for i=1:9
    subplot(3,3,i)
    optM = a2aB{i};
    a2a = [a2ad{i}(~isnan(a2aL{i})),a2aL{i}(~isnan(a2aL{i}))]';
    a2a = a2a(:,a2a(2,:)~=0);
    Plot2DHist(a2a(1,:)',a2a(2,:)',optM(1),optM(2),xlimit);
    title(['A2A: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end

figure;
xlimit = 4000;
for i=1:9
    subplot(3,3,i)
    optM = a2nB{i};
    a2n = [a2nd{i}(~isnan(a2nL{i})),a2nL{i}(~isnan(a2nL{i}))]';
    a2n = a2n(:,a2n(2,:)~=0);
    Plot2DHist(a2n(1,:)',a2n(2,:)',optM(1),optM(2),xlimit);
    title(['A2N: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end

figure;
for i=1:9
    subplot(3,3,i)
    optM = n2nB{1};
    n2n = [n2nd{i}(~isnan(n2nL{i})),n2nL{i}(~isnan(n2nL{i}))]';
    n2n = n2n(:,n2n(2,:)~=0);
    Plot2DHist(n2n(1,:)',n2n(2,:)',round(optM(1)/3),optM(2)*2,xlimit);
    title(['N2N: Culture: ', num2str(DataSet{i}.culture),' Channel: ',num2str(DataSet{i}.channel)]);
end
