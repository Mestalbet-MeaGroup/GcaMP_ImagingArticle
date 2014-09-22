function [orderedCmat,ax,nx]=DendrogramOrderMaxCosSim(MaxCosSim,ic)

a2a=MaxCosSim(size(ic,2)+1:end,size(ic,2)+1:end);
n2n=MaxCosSim(1:size(ic,2),1:size(ic,2));
a2n=MaxCosSim(1:size(ic,2),size(ic,2)+1:end);

[Na2a,ax,~]=DendrogramOrderMatrix2(a2a+a2a'+eye(size(a2a)));
[Nn2n,nx,~]=DendrogramOrderMatrix2(n2n+n2n'+eye(size(n2n)));
orderedCmat = zeros(size(MaxCosSim));
orderedCmat(1:numel(nx),1:numel(nx))=Nn2n;
orderedCmat(numel(nx)+1:numel(nx)+numel(ax),numel(nx)+1:numel(nx)+numel(ax))=Na2a;
orderedCmat(1:numel(nx),numel(nx)+1:numel(nx)+numel(ax))=a2n(nx,ax);
orderedCmat(numel(nx)+1:numel(nx)+numel(ax),1:numel(nx))=a2n(nx,ax)';
end