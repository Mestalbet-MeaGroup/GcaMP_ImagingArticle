function loadcult(which)
% global t ic mask traces time MaxCosSim
load('CorrDistributions.mat', 'MaxCosSim')
load('DataSet_GFAP_GcAMP6_withSchematic_withMask_withLags_ParCor_FullSet2.mat', 'DataSet')
t= DataSet{which}.t;
ic= DataSet{which}.ic;
mask = DataSet{which}.mask;
traces = DataSet{which}.dfTraces;
time=DataSet{which}.dfTime;
MaxCosSim = MaxCosSim{which};
a2a=MaxCosSim(size(ic,2)+1:end,size(ic,2)+1:end);
n2n=MaxCosSim(1:size(ic,2),1:size(ic,2));
assignin('base', 't', t);
assignin('base', 'ic', ic);
assignin('base', 'mask', mask);
assignin('base', 'traces', traces);
assignin('base', 'time', time);
assignin('base', 'MaxCosSim', MaxCosSim);
assignin('base', 'a2a', a2a);
assignin('base', 'n2n', n2n);
end