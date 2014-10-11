%---Figure 2---%
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2_ManSBs_withTrim.mat')
%Definitions:
%---Latency: Time between intercept with mean line and next intercept
%---Displacement: Time between burst start and astro peak
%---Amplitude: df/f value of trace
%% Subplot 1: Peaks near bursts by amplitude
k=6;
[PeakTypeBurst,PeakTypeOther,on,nc,DispNear,DispFar]=ClassifyAstroPeaks(DataSet{k}.Trim.t,DataSet{k}.Trim.ic,DataSet{k}.dfTraces,DataSet{k}.dfTime,[DataSet{k}.Trim.bs,DataSet{k}.sbs]);
amps  = nanmax(PeakTypeBurst,[],1);
numlevs=3;
ranks = otsu(amps,numlevs);
for i=1:numlevs
    subplot(3,1,i)
    plot(PeakTypeBurst(:,ranks==i),'.-');
    ylim([0,max(amps)]);
    xlim([0,100]);
end
    
%% Subplot 2: Peaks near bursts by latency
PeakTypeBurst = PeakTypeBurst(:,ranks>1);
time = [0:size(PeakTypeBurst,1)-1]./DataSet{k}.fs;
% range = 40:90;
% for j = 1:size(PeakTypeBurst,2)
%     latency(:,j) = trapz(time(range),zscore(PeakTypeBurst(range,j))+abs(min(zscore(PeakTypeBurst(range,j)))));
% end
for j = 1:size(PeakTypeBurst,2)
    temp = PeakTypeBurst(:,j)./max(PeakTypeBurst(:,j));
    newpk(:,j)=temp;
    start = find(temp>mean(temp),1,'First');
    stop = find(temp(start:end)<mean(temp),1,'First')+start;
    range = start:stop;
    if ~isempty(range)
        %         latency(j) = trapz(time(range),temp(range)+abs(min(temp(range))));
        latency(j) = numel(range);
    else
        latency(j)=nan;
    end
end
numlevs=3;
ranks = otsu(latency,numlevs);
for i=1:numlevs
    subplot(3,1,i)
    plot(zscore(PeakTypeBurst(:,ranks==i)),'.-');
    ylim([-1,5]);
    xlim([0,100]);
end

%% Subplot 3: Amplitude by Displacement from burst
% [pi,~,~,displacement]=FindPeaksNearestBursts(DataSet{k}.dfTraces,DataSet{k}.dfTime,DataSet{k}.bs);
peaks = cat(2,PeakTypeBurst,PeakTypeOther);
amplitudes=max(peaks,[],1);
displacement= cat(2,cell2mat(DispNear),cell2mat(DispFar))./DataSet{k}.fs;
scatter(displacement,amplitudes);
xlabel('Displacement [sec]');
ylabel('Amplitudes [a.u.]');

%% Subplot 4: Latency by Displacement from burst
for j = 1:size(peaks,2)
    temp = peaks(:,j)./max(peaks(:,j));
    start = find(temp>mean(temp),1,'First');
    stop = find(temp(start:end)<mean(temp),1,'First')+start;
    range = start:stop;
    if ~isempty(range)
        latency(j) = numel(range)/DataSet{k}.fs;
    else
        latency(j)=nan;
    end
end
scatter(displacement,latency);
xlabel('Displacement [sec]');
ylabel('Latency [sec]');

%% Subplot 5: Latency vs. Amplitude
scatter(latency,amplitudes);
xlabel('Latency [sec]');
ylabel('Amplitudes [a.u.]');
hist3([latency;amplitudes]',nsOPTBINS([latency;amplitudes]));

%% Subplot 6: Distribution of Amplitudes
dispHIST(amplitudes, OPTBINS(amplitudes,250));
%% Subplot 7: Distribution of Latencies
figure; 
dispHIST(latency, OPTBINS(latency,250));
%% Subplot 8: Distribution of Displacements
dispHIST(displacement, OPTBINS(displacement,250));
