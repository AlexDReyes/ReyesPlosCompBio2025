%Scripts for re-creating Figs 1-5 of manuscript
% Simulations for sustained stimulus; transient stimulus are in another script

%Each script recreates figures in manuscript; data saved
%common functions near bottom of page.

%**NOTE: highlight segment you wish to run and push "Run Selection" in Editor tab;  if you push "Run",  it will run the entire program**

%% ---------Fig. 2B of manuscript---------
totalTime=200;delT=0.01;Rec_length=round(totalTime/delT);time=delT*(1:Rec_length); %time related variables in ms
nE=250;%number of excitatory inputs
kn=0.2; %ratio of inhibitory to excitatory neurons; must be between 0-1
gEmax=0.0103; gImax=gEmax; %synaptic current peaks (current clamp); 250µV PSPs
kq=gImax/gEmax;%ratio of I to E charge or conductance; must be between 0-1;
rE=50; %rate of afferent input (Poisson process) to reference cell in Hz; 
pE=0.85; %excitatory input probability to reference cell
pI=kn*kq; %effective inhibitory probability to reference cell
pEtoI=0.35; %probability of inputs to I cells; controls input ratef
pEavg=zeros(1,Rec_length);pIavg=pEavg;pNavg=pEavg;
spkE=zeros(1,Rec_length);histE=spkE;


%-UNCOMMENT DESIRED FIGURE----------
% %----Fig 2Bi--------
% pE=0.35;
% piEtoI=0.35;
% kn=0.3;
%gMode=1; % current based, switch to 2 for conductance
% sweeps=1000;
% 
% %----Fig 2Bii--------
% pE=0.55;
% piEtoI=0.35;
% kn=0.4;
%gMode=1; % current based, switch to 2 for conductance
% sweeps=1000;

%----Fig 2Biii--------
pE=0.85;
piEtoI=0.35;
kn=0.68;
gMode=1; % current based, switch to 2 for conductance
sweeps=1000;

%------synaptic parameters
if gMode==1  %current based
    gEmax=0.0103; gImax=gEmax; %synaptic current peaks (current clamp); 250µV PSPs
end
if gMode==2  %conductance based
    gEmax=1.47e-4; %for 250µV EPSP
    gImax=10.45e-4; %for -250µV IPSP
end

tic
for sw=1:sweeps
    [pEt,pIt,rI,vE]=go(kn,kq,pE,piEtoI,rE,nE,gMode,gEmax,gImax,delT,Rec_length);
    pEavg=pEavg+pEt;
    pIavg=pIavg+pIt;
    pNavg=pNavg+(pEt.*(1-pIt));
    spkE((vE>=-10))=1;
    histE=histE+spkE;
end
pEavg=pEavg/sweeps;pIavg=pIavg/sweeps;pNavg=pNavg/sweeps;
toc
figure(1);clf;hold on;
subplot(2,1,1);plot(time,histE);xlabel('time (ms)','FontSize',14);ylabel('counts','FontSize',14)
subplot(2,1,2);hold on;plot(time,pEavg,'b');plot(time,pIavg,'r');plot(time,pNavg,'k');xlabel('time (ms)','FontSize',14);ylabel ('p','FontSize',14);lgd=legend('pE','pI','pNet');lgd.Title.FontSize = 14;

%%    Figure 3B
gMode=1; %1 is current clamp; 2 is conductance
totalTime=1000;delT=0.01;Rec_length=round(totalTime/delT);time=delT*(1:Rec_length); %time related variables in ms
%-----------constants----------------------
arStart=round(0.2*Rec_length);arEnd=round(0.8*Rec_length);

%----set up parameters to vary--------
sweeps=100;
nTrialsE=20;
gEmax=0.0103; gImax=gEmax; %synaptic current peaks (current clamp); 250µV PSPs
probE=linspace(0.1,1,nTrialsE);
piEtoI=[0 0.4 0.5 0.6 0.7];
nTrialsI=numel(piEtoI);
k=0.2;kq=1;


rE=50*ones(1,Rec_length); %firing rate of afferents in Hz
ntheta=60;
dM=zeros(nTrialsI,nTrialsE,6); %data Matrix
tic
figure(300);clf;hold on;
for jj=1:nTrialsI
    for j=1:nTrialsE
        [prbE,prbEstd,meanE,stdE,prbI,prbIstd,meanI,stdI,meanF,stdF]=goB(k,kq,probE(j),piEtoI(jj),rE(j),gEmax,gImax,gMode,sweeps,arStart,arEnd);
         ['nTrialsI: ' num2str(jj) ' nTrialsE: ' num2str(j) ' probE: ' num2str(prbE) ' probI: ' num2str(prbI) ' freqE: ' num2str(meanF) ]
        dM(jj,j,1)=prbE;
        dM(jj,j,2) = prbI; 
        dM(jj,j,3) = meanE;
        dM(jj,j,4)=meanI;
        dM(jj,j,5)=meanF;
        dM(jj,j,6)=stdF;
     %  ['pE is:' num2str(probE(j)) ' pI is: ' num2str(Iavg) ' piEtoI: ' num2str(piEtoI(jj))]
    end
    errorbar(dM(jj,:,1),dM(jj,:,5),dM(jj,:,6),'-sqk','MarkerFace','k');axis([0 1 0 200]);xlabel('pE','FontSize',20);ylabel('firing rate','FontSize',20);
end
save('Fig3B.mat','dM');
toc


%%     Figure 3C
gMode=1;%1 is current clamp; 2 is conductance
totalTime=1000;delT=0.01;Rec_length=round(totalTime/delT);time=delT*(1:Rec_length); %time related variables in ms
arStart=round(0.2*Rec_length);arEnd=round(0.8*Rec_length);
%-----------constants----------------------
gEmax=0.0103; gImax=gEmax; %synaptic current peaks (current clamp); 250µV PSPs


%----set up parameters to vary--------
sweeps=100;
nTrialsE=20;
probE=linspace(0.0,1,nTrialsE);
kScalar=[0 0.2 0.4 0.6 0.8];
nTrialsK=numel(kScalar);
rE=50;
k=1;kq=1;
piEtoI=0.35;

dM=zeros(nTrialsK,nTrialsE,6);
tic
 figure(300);clf;hold on;
for jj=1:nTrialsK
    k=kScalar(jj)*probE;
    for j=1:nTrialsE
        [prbE,prbEstd,meanE,stdE,prbI,prbIstd,meanI,stdI,meanF,stdF]=goB(k(j),kq,probE(j),piEtoI,rE,gEmax,gImax,gMode,sweeps,arStart,arEnd);
         ['nTrialsI: ' num2str(jj) ' nTrialsE: ' num2str(j) ' probE: ' num2str(prbE) ' probI: ' num2str(prbI) ' k: ' num2str(k(j)) ' freqE: ' num2str(meanF) ]
        dM(jj,j,1)=prbE;
        dM(jj,j,2) = prbI; 
        dM(jj,j,3) = meanE;
        dM(jj,j,4)=meanI;
        dM(jj,j,5)=meanF;
        dM(jj,j,6)=stdF;
    end
     errorbar(dM(jj,:,1),dM(jj,:,5),dM(jj,:,6),'-sqk','MarkerFace','k');axis([0 1 0 200]);xlabel('pE','FontSize',20);ylabel('firing rate','FontSize',20)
end
save('Fig3C.mat','dM');
toc
%% Figure 3D
gMode=1;
totalTime=1000;delT=0.01;Rec_length=round(totalTime/delT);time=delT*(1:Rec_length); %time related variables in ms
arStart=round(0.2*Rec_length);arEnd=round(0.8*Rec_length);
%-----------constants----------------------


%----set up parameters to vary--------
sweeps=100;
nTrialsE=20;
probE=linspace(0.0,1,nTrialsE);
pEtoI=[0 0.4 0.5 0.6 0.7 0.8];
nTrialspEtoI=numel(pEtoI);
gEmax=0.0103; gImax=gEmax; %synaptic current peaks (current clamp); 250µV PSPs
k=0.2;kq=1;
rE=50;

dM=zeros(nTrialspEtoI,nTrialsE,6);

tic
 figure(300);clf;hold on;
for jj=1:nTrialspEtoI
    piEtoI=pEtoI(jj)*probE;
    for j=1:nTrialsE
        [prbE,prbEstd,meanE,stdE,prbI,prbIstd,meanI,stdI,meanF,stdF]=goB(k,kq,probE(j),piEtoI(j),rE,gEmax,gImax,gMode,sweeps,arStart,arEnd);
        ['nTrialspEtoI: ' num2str(jj) ' nTrialsE: ' num2str(j) ' probE: ' num2str(prbE) ' probI: ' num2str(prbI) ' freqE: ' num2str(meanF) ]
        dM(jj,j,1)=prbE;
        dM(jj,j,2) = prbI; 
        dM(jj,j,3) = meanE;
        dM(jj,j,4)=meanI;
        dM(jj,j,5)=meanF;
        dM(jj,j,6)=stdF;
    end
    pause(0.1);
     errorbar(dM(jj,:,1),dM(jj,:,5),dM(jj,:,6),'-sqk','MarkerFace','k');axis([0 1 0 200]);xlabel('pE','FontSize',20);ylabel('firing rate','FontSize',20)
end
    save('Fig3D.mat','dM');
toc

%% ----Figure 4B-------------
gMode=1; %1 is current clamp; 2 is conductance
totalTime=1000;delT=0.01;Rec_length=round(totalTime/delT);time=delT*(1:Rec_length); %time related variables in ms
arStart=round(0.2*Rec_length);arEnd=round(0.8*Rec_length);
%-----------constants----------------------
gEmax=0.0103; gImax=gEmax; %synaptic current peaks (current clamp); 250µV PSPs


%----set up parameters to vary--------
sweeps=100;
nTrialsE=50;
x=linspace(2,8,nTrialsE);
xctr=round(mean(x));sigX=1.5;
k=0.2;kq=1;
pEtoI=[0 0.275 0.3]; %synaptic efficacy from I cells to reference cell
nTrialspEtoI=numel(pEtoI);
aScalar=0.35;
rE=50;

dM=zeros(nTrialspEtoI,nTrialsE,6);

tic
figure(300);clf;hold on;
probE=aScalar*exp(-0.5*((x-xctr)/sigX).^2);
for jj=1:nTrialspEtoI
    for j=1:nTrialsE
        [prbE,prbEstd,meanE,stdE,prbI,prbIstd,meanI,stdI,meanF,stdF]=goB(k,kq,probE(j),pEtoI(jj),rE,gEmax,gImax,gMode,sweeps,arStart,arEnd);
        ['nTrialspEtoI: ' num2str(jj) ' nTrialsE: ' num2str(j) ' probE: ' num2str(prbE) ' probI: ' num2str(prbI) ' freqE: ' num2str(meanF) ]
        dM(jj,j,1)=prbE;
        dM(jj,j,2) = prbI; 
        dM(jj,j,3) = meanE;
        dM(jj,j,4)=meanI;
        dM(jj,j,5)=meanF;
        dM(jj,j,6)=stdF;
    end
     figure(300);errorbar(x,dM(jj,:,5),dM(jj,:,6),'-sqk','MarkerFace','k');axis([min(x) max(x) -2 60]);xlabel('sensory Feature','FontSize',20);ylabel('firing rate','FontSize',20)
end
save('Fig4B.mat','dM');
toc

%%   Figure 4C
gMode=1;
totalTime=1000;delT=0.01;Rec_length=round(totalTime/delT);time=delT*(1:Rec_length); %time related variables in ms
arStart=round(0.2*Rec_length);arEnd=round(0.8*Rec_length);
%-----------constants----------------------
gEmax=0.0103; gImax=gEmax; %synaptic current peaks (current clamp); 250µV PSPs

%----set up parameters to vary--------
sweeps=100;
nTrialsE=20;
x=linspace(2,8,nTrialsE);
xctr=round(mean(x));sigX=1.5;
kScalar=[0 0.2 0.4];
nTrialsK=numel(kScalar);
aScalar=0.35;
dM=zeros(nTrialsK,nTrialsE,6);
rE=50;
k=1;
piEtoI=0.35;
tic
 figure(300);clf;hold on;
for jj=1:nTrialsK
    probE=aScalar*exp(-0.5*((x-xctr)/sigX).^2);
    k=kScalar(jj)*probE;
    for j=1:nTrialsE
        [prbE,prbEstd,meanE,stdE,prbI,prbIstd,meanI,stdI,meanF,stdF]=goB(k(j),kq,probE(j),piEtoI,rE,gEmax,gImax,gMode,sweeps,arStart,arEnd);
         ['nTrialsI: ' num2str(jj) ' nTrialsE: ' num2str(j) ' probE: ' num2str(prbE) ' probI: ' num2str(prbI) ' k: ' num2str(k(j)) ' freqE: ' num2str(meanF) ]
        dM(jj,j,1)=prbE;
        dM(jj,j,2) = prbI; 
        dM(jj,j,3) = meanE;
        dM(jj,j,4)=meanI;
        dM(jj,j,5)=meanF;
        dM(jj,j,6)=stdF;
    end
     figure(300);errorbar(x,dM(jj,:,5),dM(jj,:,6),'-sqk','MarkerFace','k');axis([min(x) max(x) -2 60]);xlabel('sensory Feature','FontSize',20);ylabel('firing rate','FontSize',20)
end
save('Fig4C.mat','dM');
toc

%% Figure 4D
gMode=1;
totalTime=1000;delT=0.01;Rec_length=round(totalTime/delT);time=delT*(1:Rec_length); %time related variables in ms
arStart=round(0.2*Rec_length);arEnd=round(0.8*Rec_length);
%-----------constants----------------------
gEmax=0.0103; gImax=gEmax; %synaptic current peaks (current clamp); 250µV PSPs

%----set up parameters to vary--------
sweeps=100;
nTrialsE=20;
x=linspace(2,8,nTrialsE);
xctr=round(mean(x));sigX=1.5;
pEtoI=[0.7 0.8 0.9];
nTrialspEtoI=numel(pEtoI);
k=0.2;kq=1;
rE=50;
aScalar=0.35;

dM=zeros(nTrialspEtoI,nTrialsE,6);

tic
 figure(300);clf;hold on;
 probE=aScalar*exp(-0.5*((x-xctr)/sigX).^2);
for jj=1:nTrialspEtoI
    piEtoI=pEtoI(jj)*probE;
    for j=1:nTrialsE
        [prbE,prbEstd,meanE,stdE,prbI,prbIstd,meanI,stdI,meanF,stdF]=goB(k,kq,probE(j),piEtoI(j),rE,gEmax,gImax,gMode,sweeps,arStart,arEnd);
        ['nTrialspEtoI: ' num2str(jj) ' nTrialsE: ' num2str(j) ' probE: ' num2str(prbE) ' probI: ' num2str(prbI) ' freqE: ' num2str(meanF) ]
        dM(jj,j,1)=prbE;
        dM(jj,j,2) = prbI; 
        dM(jj,j,3) = meanE;
        dM(jj,j,4)=meanI;
        dM(jj,j,5)=meanF;
        dM(jj,j,6)=stdF;
    end
    pause(0.1);
     errorbar(x,dM(jj,:,5),dM(jj,:,6),'-sqk','MarkerFace','k');axis([min(x) max(x) -2 60]);xlabel('sensory Feature','FontSize',20);ylabel('firing rate','FontSize',20)
end
    save('Fig4D.mat','dM');
toc


  %%   Figure 5
totalTime=120;delT=0.01;Rec_length=round(totalTime/delT);time=delT*(1:Rec_length); %time related variables in ms
histE=zeros(1,Rec_length);
rE=50;rI=rE;gEmax=0.0103; gImax=gEmax; 

%----parameters to vary
sweeps=1000;
pE=1;% steady-state excitatory probability
delE=5;delI=5; % relative delay (onset);example: if delE=5, delI=7 then inhibition occurs 2 ms after excitation
rmpE=round(20);rmpI=rmpE; %ramp duration in ms
nE=250;nI=nE;

[pEt,pIt]=ramp(delE,delI,rmpE,rmpI,totalTime);
pEt=pE*smooth(pEt,500);pEt=pEt';
pIt=pE*smooth(pIt,500);pIt=pIt';
pNt=pEt.*(1-pIt);

tic
for sw=1:sweeps
   gE=gEmax*genPSPtrain(pEt,rE,nE,totalTime,delT);
   gI=gImax*genPSPtrain(pEt.*pIt,rI,nI,totalTime,delT);
   v=LIF(gE,gI,1,Rec_length,delT,1);
   histE(find(v>=0))= histE(find(v>=0))+1;
end
figure(1);clf;hold on;subplot(2,1,1);hold on;plot(time,pEt,'b','lineWidth',2);plot(time,pIt,'--r','lineWidth',3);plot(time,pNt,'k','lineWidth',2);axis([0 120 0 1]);
    xlabel('time (ms)','FontSize',20);ylabel('Prob','FontSize',20);lgd=legend('pE','pI','pNet');lgd.Title.FontSize = 14;
    subplot(2,1,2);stem(time,histE,'k','Marker','none');xlabel('time (ms)','FontSize',20);ylabel('counts','FontSize',20);xlim([0 120]);
h=figure(1);savefig(h,'fig5');
  toc      


%%  common functions
function [gdIn]=genPSPtrain(probIn,r,n,totalTime,delT)
Rec_length=round(totalTime/delT);
        ptau=2.0;impLength=round(25/delT);t=(1:impLength);
        uPSP=(delT*t/ptau).*exp(1-t*delT/ptau);
        gdIn=zeros(1,numel(probIn));gdI=gdIn;gImpE=zeros(1,round(25/delT));

        %train=binornd(nE,tR*delT*0.001*probIn);
        train=poissrnd(n*r*delT*0.001*probIn); 

        gdIn = conv(train,uPSP);
        gdIn(Rec_length+1:numel(gdIn))=[];
end

function [envE,envI]=ramp(delE,delI,rmpE,rmpI,totalTime)
    delT=0.01;rampDur=rmpE;del=round(delE/delT);
    Rec_length=round(totalTime/delT);envE=zeros(1,Rec_length);
    dur=round((totalTime-2*rampDur-2*delE)/delT);rLength=round(rampDur/delT);m=1/rLength;xEnd = rLength+del+1+dur-100;
    envE(del:rLength+(del))=(1/rLength)*(del:del+rLength)-(1/rLength)*del;
    envE(rLength+del+1:xEnd)=1;
    envE(xEnd:rLength+xEnd-1)=-m*(xEnd:rLength+xEnd-1)+m*(xEnd+rLength);

    delT=0.01;rampDur=rmpI;del=round(delI/delT);diffDel=(round(delE/delT)-round(delI/delT));
    Rec_length=round(totalTime/delT);envI=zeros(1,Rec_length);
    dur=round((totalTime-2*rampDur-2*delI)/delT)-diffDel;rLength=round(rampDur/delT);m=1/rLength;
    xEnd = rLength+del+1+dur-100-diffDel;
    envI(del:rLength+(del))=(1/rLength)*(del:del+rLength)-(1/rLength)*del;
    envI(rLength+del+1:xEnd)=1;
    envI(xEnd:rLength+xEnd-1)=-m*(xEnd:rLength+xEnd-1)+m*(xEnd+rLength);

end

function [pEt,pIt,rI,vE]=go(kn,kq,pE,pEtoI,rE,nE,gMode,gEmax,gImax,delT,Rec_length)
  %--constants,vectors, matrices-------       
     nI=kn*nE;  %number of inhibitory cells
     vE=zeros(1,Rec_length); vI=vE;%voltage traces
     spkI=zeros(1,Rec_length);
     totalTime=delT*Rec_length;
     %---synaptic parameters-------
     gE=zeros(1,Rec_length);gI=kq*gE; %synaptic conductances 
     pEt=zeros(1,Rec_length); pIt=pEt;%probability traces
     %-----make template synaptic current; alpha function------
     ptau=2;t=(1:Rec_length); uPSP=1*(delT*t/ptau).*exp(1-t*delT/ptau);

     %------make a bank of inhibitory inputs to E cells---
     %simulate the activities of I cells and store in matrices
        numWaves=3*nI; %3x as many as needed; goes slower if increases
        spkIbank=zeros(numWaves,Rec_length); %matrices to store spikes
        for nn=1:numWaves 
            train=poissrnd(pEtoI*rE*nE*delT*0.001,[1 Rec_length]); %when pE=1
            gE=gEmax*conv(train,uPSP);gE(Rec_length+1:numel(gE))=[];
            vI=LIF(gE,0*gI,gMode,Rec_length,delT,1); %with (1) or without (0) spikes
            spkIbank(nn,(find(vI>=0)))=spkIbank(nn,(find(vI>=0)))+1;
        end
rI=1000*sum(sum(spkIbank(:,:)))/(numWaves*totalTime); %average firing rate of inhibitory cells
  %-------calculate E and I inputs to reference cell-----------
    gI=gI*0;
    lambda = rE*nE*delT*0.001;
    train=poissrnd(pE*lambda,[1 Rec_length]);
    cPSP=conv(train,uPSP);cPSP(Rec_length+1:numel(cPSP))=[];
    pEt=cPSP/(sum(uPSP)*lambda);% calculate excitatory probability time course
    gE=gEmax*cPSP;
    Ilist=randperm(numWaves,nI); %pick random set of nI inhibitory inputs from bank
    for ii=1:numel(Ilist)
        spkI=spkI+spkIbank(Ilist(ii),:);
    end
    cPSP=conv(spkI,uPSP);cPSP(Rec_length+1:numel(cPSP))=[];
    pIt=cPSP/(sum(uPSP)*lambda);
    gI=gImax*cPSP;%conv(spkI,kq*uPSP);gI(Rec_length+1:numel(gI))=[];
     vE=LIF(gE,pE*gI,gMode,Rec_length,delT,1);
     %nEspks=sum(vE>-10); %count number of evoked spikes
end

function [prbE,prbEstd,meanE,stdE,prbI,prbIstd,meanI,stdI,meanF,stdF]=goB(kval,kq,prbE,pIX, tR,gEmax,gImax,gMode,sweeps,arStart,arEnd)
        nE=250;nI=nE;delT=0.01;totalTime=1000;Rec_length=round(totalTime/delT); IRh=0.2;q=0.0557;tau=10;
        lambda=nE*tR*delT*0.001; %used for Poisson train
        vE=zeros(1,Rec_length);vI=vE;freqN=zeros(1,Rec_length);
        gE=zeros(1,Rec_length);gI=gE;
        ptau=2;t=(1:Rec_length); uPSP=(delT*t/ptau).*exp(1-t*delT/ptau); imp=exp(-t/500);imp=imp/sum(imp);
        knI=round(kval*nI);
        kswps=3*knI; spkIbank=zeros(kswps,Rec_length);   
        nIspks=0;nMaxspks=0;
        for nn=1:kswps %rI
        train=poissrnd(nE*tR*delT*0.001*pIX,1,Rec_length);
        gE = gEmax*conv(train,uPSP);gE(Rec_length+1:numel(gE))=[];
            vI=LIF(gE,0*gI,1,Rec_length,delT,1) ;
            spkIbank(nn,(find(vI>=0)))=spkIbank(nn,(find(vI>=0)))+1;
            nMaxspks=nMaxspks+sum(vI(arStart:arEnd)>-10);
            if(sum(vI>-10))>0
                nIspks=nIspks+1;
            end
        end
        nMaxspks=nMaxspks/kswps;
        rI=nMaxspks/((arEnd-arStart)*delT/1000);
        ppI=nIspks/kswps;
        %-----make unitary PSP template--------
        ptau=2;t=(1:Rec_length);
        uPSP=(delT*t/ptau).*exp(1-t*delT/ptau);
        nEspks=0;nEspksPred=0;freqIavg=0;spkMatrix=zeros(4,sweeps);
        q=sum(gEmax*uPSP)*delT;
        fspks=zeros(1,sweeps);
        for sw=1:sweeps
            gIhisto=zeros(1,Rec_length);
            Ilist=randperm(kswps,knI);
            for nn=1:knI
                gIhisto=gIhisto+spkIbank(Ilist(nn),:);
            end
            cPSP=conv(gIhisto,kq*uPSP);cPSP(Rec_length+1:numel(cPSP))=[];
            pIt=cPSP/(sum(uPSP)*lambda);
            gI=gImax*cPSP;%conv(spkI,kq*uPSP);gI(Rec_length+1:numel(gI))=[]; 
            spkMatrix(3,sw)=mean(pIt(arStart:arEnd));  %mean I prob
            spkMatrix(4,sw)=mean(gI(arStart:arEnd)); %mean I current
            gI=gI*prbE;%conditioned on pE
           
            train=poissrnd(lambda*prbE,1,Rec_length);
            cPSP=conv(train,uPSP);cPSP(Rec_length+1:numel(cPSP))=[];
            pEt=cPSP/(sum(uPSP)*lambda);% calculate excitatory probability time course
            gE = gEmax*cPSP;%conv(train,uPSP);gE(Rec_length+1:numel(gE))=[];
            spkMatrix(1,sw)=mean(pEt(arStart:arEnd)); %mean E prob
            spkMatrix(2,sw)=mean(gE(arStart:arEnd));  %mean E current
            vE=LIF(gE,gI,1,Rec_length,delT,1);%figure(101);clf;hold on;plot(vE,'b');
            nEspks=sum(vE>-10);
            spkMatrix(5,sw)=1000*nEspks/totalTime;  %firing rate in Hz
        end
        prbE=mean(spkMatrix(1,:));prbEstd=std(spkMatrix(1,:));
        meanE=mean(spkMatrix(2,:))/sweeps;stdE=std(spkMatrix(2,:));
        prbI=mean(spkMatrix(3,:));prbIstd=std(spkMatrix(3,:));
        meanI=mean(spkMatrix(4,:))/sweeps;stdI=std(spkMatrix(4,:));
        meanF=mean(spkMatrix(5,:));stdF=std(spkMatrix(5,:));
end

function v=LIF(gexc,ginh,gMode,Rec_length,delT,spkOn) 
    tau=10;R=75;El=-70;vTh=-55; %LIF parameters
    v=zeros(1,Rec_length)+El;Isyn=0*v;
   i=2;
    while (i<=Rec_length)
        if gMode==2   %conductance mode
             Isyn(i)=gexc(i)*(v(i-1)-0)+ginh(i)*(v(i-1)+80);
        end
        if gMode==1   %current clamp mode
            Isyn(i)=-(gexc(i)-ginh(i));
        end
        delV=(-(v(i-1)-El)-R*(Isyn(i)))*(delT/tau);  
        v(i)=v(i-1)+delV;    
        if spkOn==1  %generate spikes?
            if (v(i)>vTh)
                v(i)=0;
                v(i+1)=El;
                i=i+1;
            end
        end
        i=i+1;
    end 
end
