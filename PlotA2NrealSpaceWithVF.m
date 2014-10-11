% plot a2n links on realspace including n2n of VF elecs 

howmany = 500;
[e1,e2,eval,a1,a_e2,aval,ch,mask]=ParseConnectionValuesWithVF(6,howmany);
e1=e1(howmany+1:end);
e2=e2(howmany+1:end);
eval = eval(howmany+1:end);
ind = find(eval<min(aval));
e1(ind)=[];
e2(ind)=[];
eval(ind)=[];

figure('Renderer','opengl');
PlotResultsOnMEAwithVF(e1,e2,eval,mask,ch,a1,a_e2,aval);


