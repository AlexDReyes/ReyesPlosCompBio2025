% Scripts for recreating Figure 6 in manuscript.  For transient inputs; Sustained inputs in another script
%each script for each figure

%**NOTE: highlight segment you wish to run and push "Run Selection" in Editor tab;  if you push "Run",  it will run the entire program**

%% Figure 6A
%Program for varying times of transient PSPs
gMode=0; %0 for current clamp, 1 for conductance clamp
delT=0.01;totalTime=50;Rec_length=round(totalTime/delT);ms=delT*(1:Rec_length);

if (gMode==0)
    gEmax=0.25*0.041; gImax=gEmax; %current clamp
    nFire=68;
end
if gMode==1
    gEmax=1.47*1e-4; %for 250µV EPSP
    gImax=10.45*1e-4; %for -250µV IPSP
    nFire=68;
end

%--------parameters to vary--------
sweeps=100;
nE=100;
%uncomment for desired figure
%pE=0.95;pIX=0.95;k=0.2;Edelay=round(10/delT);Esig=round(1/delT);Idelay=round(10/delT);Isig=round(1/delT);%Figure 6Ai
pE=0.75;pIX=0.75;k=0.2;Edelay=round(10/delT);Esig=round(1/delT);Idelay=round(10/delT);Isig=round(1/delT);%Figure 6Aii

[histE,histI,histX,vEPSP,probE,probI,probNet,fN,fMatrix]=go(pE,pIX,k,Edelay,Esig,Idelay,Isig,totalTime,delT,sweeps,nE,gMode,nFire);

probNet(find(probNet>=1))=1;probNet(find(probNet<=0))=0;
pTheta=1-binocdf(nFire,nE,probNet);

figure(200);clf;hold on;
subplot(3,1,1);hold on;stem(ms,delT*histX,'b','Marker','none');stem(ms,delT*histI,'r','Marker','none');
subplot(3,1,3);hold on;plot(ms,probE,'b','lineWidth',2);plot(ms,probI,'r','lineWidth',2);plot(ms,probNet,'k','lineWidth',2);plot(ms,pTheta,'g','lineWidth',2);ylim([0 1]);
subplot(3,1,2);hold on;plot(ms,fN,'m','lineWidth',1);bar(ms,histE,'FaceColor',[0.6 0.6 0.6]);
figure(201);clf;hold on;bar(ms,histE,'FaceColor',[0.6 0.6 0.6]);plot(ms,fN,'m','lineWidth',2);


%% Figure 6Bi
delT=0.01;totalTime=30;Rec_length=round(totalTime/delT);t=(1:Rec_length);
pEt=zeros(1,Rec_length);hD=pEt;
gEmax=0.25*0.041; gImax=gEmax; %unitary PSC amplitude for 250µV PSP
nFire=75; %mininum number of inputs to evoke firing

%----parameters to vary
sweeps=1000;
nStim=20;%number of points
nE=100;
piEtoI=[0 0.7 0.8 0.9 1]; %synaptic efficacy values
pStim=linspace(0.4,1,nStim);nProbs=numel(pStim); %range of pE values
k=0.2; %ratio of inhib to exc
Edelay=round(10/delT);Esig=round(1/delT);Idelay=round(10/delT);Isig=round(1/delT);%relative delays between E and I

figure(300);clf;hold on;figure(301);clf;hold on;
tic
meanHistE=zeros(numel(piEtoI),nProbs);stdHistE=meanHistE;meanHistN=meanHistE;stdHistN=meanHistE;meanHistI=meanHistE;stdHistE=meanHistE;stdHistI=meanHistE;
for ii=1:numel(piEtoI)
    for jj=1:nProbs
       [histE,histI,histX,vEPSP,probE,probI,probNet,fN,fMatrix]=go(pStim(jj),piEtoI(ii), k,Edelay,Esig,Idelay,Isig,totalTime,delT,sweeps,nE,0,nFire); 
        meanHistE(ii,jj)=fMatrix(1);stdHistE(ii,jj)=fMatrix(2);
        meanHistN(ii,jj)=fMatrix(3);stdHistN(ii,jj)=fMatrix(4);
        meanHistI(ii,jj)=fMatrix(5);stdHistI(ii,jj)=fMatrix(6);
    end
    figure(300);subplot(2,1,1);hold on;plot(pStim,meanHistN(ii,:),'-sqk','markerfacecolor','k');
    subplot(2,1,2);hold on;plot(pStim,meanHistE(ii,:),'-sqk','markerfacecolor','k');plot(pStim,meanHistI(ii,:),'-sqr','markerfacecolor','r');
    figure(301);plot(pStim,meanHistN(ii,:),'g','lineWidth',2);errorbar(pStim,meanHistE(ii,:),stdHistE(ii,:),'-sqk','markerfacecolor','k');
end
h(1)=figure(300);h(2)=figure(301);
save('Fig6Bi','meanHistE','stdHistE','meanHistI','stdHistI','meanHistN','stdHistN')
toc

%%  Figure 6Bii
delT=0.01;totalTime=30;Rec_length=round(totalTime/delT);t=(1:Rec_length);
gEmax=0.0103; gImax=gEmax; %unitary PSC amplitude for 250µV unitary PSP
nFire=74;

%----parameters to vary------
sweeps=1000;
nStim=20;
nE=100;
pStim=linspace(0.4,1,nStim);nProbs=numel(pStim);% number of pE points
piEtoI=0.8; %synaptic efficacy
kIval=[0 0.2 0.4 0.6 0.8 1]; %scales inhibition so increases linearly with pE
k=1; %ratio of inhib to exc
Edelay=round(10/delT);Esig=round(1/delT);Idelay=round(10/delT);Isig=round(1/delT); %relative delays between E and I

figure(300);clf;hold on;figure(301);clf;hold on;
tic
meanHistE=zeros(numel(kIval),nProbs);stdHistE=meanHistE;meanHistN=meanHistE;stdHistN=meanHistE;meanHistI=meanHistE;stdHistE=meanHistE;stdHistI=meanHistE;
 for ii=1:numel(kIval)
     for jj=1:nProbs
       [histE,histI,histX,vEPSP,probE,probI,probNet,fN,fMatrix]=go(pStim(jj),piEtoI, kIval(ii)*pStim(jj),Edelay,Esig,Idelay,Isig,totalTime,delT,sweeps,nE,0,nFire);
        meanHistE(ii,jj)=fMatrix(1);stdHistE(ii,jj)=fMatrix(2);
        meanHistN(ii,jj)=fMatrix(3);stdHistN(ii,jj)=fMatrix(4);
        meanHistI(ii,jj)=fMatrix(5);stdHistI(ii,jj)=fMatrix(6);
    end
    pause(0.01);
    figure(300);subplot(2,1,1);hold on;plot(pStim,meanHistN(ii,:),'-sqk','markerfacecolor','k');
    subplot(2,1,2);hold on;plot(pStim,meanHistE(ii,:),'-sqk','markerfacecolor','k');plot(pStim,meanHistI(ii,:),'-sqr','markerfacecolor','r');
    figure(301);hold on;plot(pStim,meanHistN(ii,:),'g','lineWidth',2);errorbar(pStim,meanHistE(ii,:),stdHistE(ii,:),'-sqk','markerfacecolor','k');
    ['ii: ' num2str(ii) ' jj: ' num2str(jj)]
 end
h(1)=figure(300);h(2)=figure(301);
save('Fig6Bii','meanHistE','stdHistE','meanHistI','stdHistI','meanHistN','stdHistN')
toc
%% Figure 6Biii

delT=0.01;totalTime=30;Rec_length=round(totalTime/delT);t=(1:Rec_length);
gEmax=0.0103; gImax=gEmax; %current clamp
nFire=73;

%--parameters to vary
sweeps=100;
nStim=20;
nE=100;
pStim=linspace(0.4,1,nStim);nProbs=numel(pStim); % number of pE points
kIval=[0 0.85 0.875 0.9 0.925 0.95 0.975 1]; %scales syn efficacy to I cell so it varies linearly with pE
k=1; %ratio of inhib to exc
Edelay=round(10/delT);Esig=round(1/delT);Idelay=round(10/delT);Isig=round(1/delT); %relative delays between E and I

figure(300);clf;hold on;figure(301);clf;hold on;
tic
meanHistE=zeros(numel(kIval),nProbs);stdHIstE=meanHistE;meanHistN=meanHistE;stdHistN=meanHistE;meanHistI=meanHistE;stdHistE=meanHistE;stdHistI=meanHistE;
 for ii=1:numel(kIval)
     for jj=1:nProbs
       [histE,histI,histX,vEPSP,probE,probI,probNet,fN,fMatrix]=go(pStim(jj),kIval(ii)*pStim(jj), 1,Edelay,Esig,Idelay,Isig,totalTime,delT,sweeps,nE,0,nFire);
        meanHistE(ii,jj)=fMatrix(1);stdHistE(ii,jj)=fMatrix(2);
        meanHistN(ii,jj)=fMatrix(3);stdHistN(ii,jj)=fMatrix(4);
        meanHistI(ii,jj)=fMatrix(5);stdHistI(ii,jj)=fMatrix(6);
    end
    pause(0.01);
    figure(300);subplot(2,1,1);hold on;plot(pStim,meanHistN(ii,:),'-sqk','markerfacecolor','k');
    subplot(2,1,2);hold on;plot(pStim,meanHistE(ii,:),'-sqk','markerfacecolor','k');plot(pStim,meanHistI(ii,:),'-sqr','markerfacecolor','r');
    figure(301);hold on;plot(pStim,meanHistN(ii,:),'g','lineWidth',2);errorbar(pStim,meanHistE(ii,:),stdHistE(ii,:),'-sqk','markerfacecolor','k');
    ['ii: ' num2str(ii) ' jj: ' num2str(jj)]
end
save('Fig6Biii','meanHistE','stdHistE','meanHistI','stdHistI','meanHistN','stdHistN')
toc

%% functions
 function [histE,histI,histX,vEPSP,probE,probI,probNet,fN,fMatrix]=go(pEval,pIX,kval,Edelay,Esig,Idelay,Isig,totalTime,delT,sweeps,nE,gMode,nFire)
        %-----------constants----------------------
        
        nI=nE;nTval=60;vTh=-55;aE=0.25;aI=0.25;
        Rec_length=round(totalTime/delT); 

        %-------vectors--------
        histIavg=zeros(1,Rec_length);
        histE=zeros(1,Rec_length);
        pEt=zeros(1,Rec_length);gDavg=pEt;
        gdIavg=zeros(1,Rec_length);probEavg=gdIavg; probIavg=gdIavg;probNetavg=gdIavg;fNavg=gdIavg;
        probI=zeros(1,Rec_length);probE=probI;probNet=probI;sfN=zeros(sweeps,3);fMatrix=zeros(1,6);
        histE=histE*0;histI=zeros(1,Rec_length);
       
        %-----------LIF parameters---------
        if (gMode==0)
            gEmax=0.0103; gImax=gEmax; %current clamp
        end
        if gMode==1
            gEmax=1.47*1e-4; %for 250µV EPSP
            gImax=10.45*1e-4; %for -250µV IPSP
          %  gImax=0.35*10.45*1e-4;% will generate same firing as current clamp
        end
        gE=zeros(1,Rec_length);gI=gE;gD=gE;El=-70;
        %-----make unitary PSP template--------
        ptau=2;t=(1:Rec_length);
        uPSP=(delT*t/ptau).*exp(1-t*delT/ptau);
         pulse=zeros(1,Rec_length);pulse(Edelay)=1;
        gEPSP=conv(gEmax*uPSP,pulse);gEPSP(Rec_length+1:numel(gEPSP))=[];
        vEPSP=LIF(gEPSP,0*gI,vTh,totalTime,delT,0,gMode)-El;
        gIPSP=conv(gImax*uPSP,pulse);gIPSP(Rec_length+1:numel(gIPSP))=[];
        vIPSP=abs(LIF(0*gEPSP,gIPSP,vTh,totalTime,delT,0,gMode)-El);
        
            %-------make a bank of inhibitory trains-----
                 %-so don't have to re-calculate inhibitory trains for each sweep
        knI=round(kval*nI);
        kswps=3*knI; spkIbank=zeros(kswps,Rec_length);histXavg=zeros(1,Rec_length);
        pHat=0;
        for nn=1:kswps 
            gD=genPreG(pIX,Idelay,Isig,gD*0,nE);
            histXavg=histXavg+gD;
            gE=gEmax*conv(gD,uPSP);gE(Rec_length+1:numel(gE))=[];
            vI=LIF(gE,0*gI,vTh,totalTime,delT,1,gMode);
            b=find(vI>-10);
            if (numel(b))>0
                spkIbank(nn,(b(1)))=spkIbank(nn,(b(1)))+1;
                pHat=pHat+1;
            end
        end
        pHat=pHat/kswps;
        gD=genPreG(1,Edelay,Esig,gD*0,nE);
        gdE=conv(gD,uPSP);gdE(Rec_length+1:numel(gdE))=[];
        gE=gEmax*gdE;
        vEPSP=LIF(gE,0*gI,vTh,totalTime,delT,0,gMode)-El;
        probE=pEval*vEPSP/(0.25*nE);
        Ilist=randperm(kswps,knI); %chooses without repeat trains from I bank
        for nn=1:knI
           histI=histI+spkIbank(Ilist(nn),:);
        end
        gdIb=conv(histI,uPSP);gdIb(Rec_length+1:numel(gdIb))=[];
        gIb=gImax*gdIb;
        vIPSP=abs(LIF(0*gEPSP,gIb,vTh,totalTime,delT,0,gMode)-El);
        probI=pHat*kval*pIX*vIPSP/(0.25*nI);
        probNet=probE-pEval*probI;
        fN=calc_fN(probNet,nE,aE,nFire);  

        %-----------Start main program-----------

        for sw=1:sweeps
            histEt=0;
            knI=round(kval*nI*pHat);
            histI=zeros(1,Rec_length);
            if knI>0
                Ilist=randperm(kswps,knI); %chooses without repeat trains from I bank
                for nn=1:knI
                   histI=histI+spkIbank(Ilist(nn),:);
                end
            end
            histIavg=histIavg+histI;
       %-----calculate E and I input to E cells-------
            gD=genPreG(pEval,Edelay,Esig,gD*0,nE);
            gDavg=gDavg+gD;
            gdE=conv(gD,uPSP);gdE(Rec_length+1:numel(gdE))=[];
            gE=gEmax*gdE;
            gdIb=conv(histI,uPSP);gdIb(Rec_length+1:numel(gdIb))=[];
            gIb=gImax*(pEval.*gdIb);

            vE=LIF(gE,gIb,vTh,totalTime,delT,1,gMode);
            
            b=find(vE>-10);
            if (numel(b))>0
                histE(b(1))=histE(b(1))+1;
                histEt=1;
            end
              vEPSP=LIF(gE,0*gI,vTh,totalTime,delT,0,gMode)-El;
              probE=vEPSP/(0.25*nE);  %uPSP amplitude
              probEavg=probEavg+probE;
              vIPSP=abs(LIF(0*gEPSP,gIb,vTh,totalTime,delT,0,gMode)-El);
                probI=vIPSP/(0.25*nI);
               probIavg=probIavg+probI;
                probNet=probE-pEval*probI;
                probNetavg=probNetavg+probNet;
                probNet(probNet<0)=0;

              fN=calc_fN(probNet,nE,aE,nFire);  
              fNavg=fNavg+fN;
              sfN(sw,1)=histEt;sfN(sw,2)=sum(fN);
               sfN(sw,3)=sum(histI)/nI;

        end
        fMatrix(1)=mean(sfN(:,1));fMatrix(2)=std(sfN(:,1));
        fMatrix(3)=mean(sfN(:,2));fMatrix(4)=std(sfN(:,2));
        fMatrix(5)=mean(sfN(:,3));fMatrix(6)=std(sfN(:,3));
        probE=probEavg/sweeps;
        probI=probIavg/sweeps;
        probNet=probNetavg/sweeps;
        fN=fNavg/sweeps;
        histX=gDavg/sweeps;
        histE=histE/sweeps;
        histI=histIavg/(sweeps);
    end

    function fN=calc_fN(probNet,nE,aE,nFire)   
    vTheta=nFire*0.25;
        vAmp = nE*aE*probNet; vSTD = aE*sqrt(nE*probNet.*(1-probNet));       
        pTheta=1-cdf('Normal',vTheta,vAmp,vSTD);
        numN=numel(pTheta);
        pB=0;fN=zeros(1,numN);
        for i=3:numN
            j=i-2;
            pB=(1-pTheta(i-1));
            while j>0
                pB=pB*(1-pTheta(j));
                j=j-1;
            end
            fN(i)=pTheta(i)*pB;
        end
    end

    function gD=genPreG(pp,Edelay,Esig,gD,nE)
        for i=1:nE
            if unifrnd(0,1)<=pp
                d=0;
               while d<=0
                    d=round(normrnd(Edelay,Esig));
               end
               gD(d)=gD(d)+1;
            end
        end
    end

    function v=LIF(gexc,ginh,thetaval,totalTime,delT,spkOn,gMode)
        Rec_length=round(totalTime/delT);
        tau=10;R=75;El=-70;
        v=zeros(1,Rec_length)+El;Isyn=v;
       i=2;
        while (i<=Rec_length)
            if gMode==1
                 Isyn(i)=gexc(i)*(v(i-1)-0)+ginh(i)*(v(i-1)+80);
            end
            if gMode==0
                Isyn(i)=-(gexc(i)-ginh(i));
            end
            delV=(-(v(i-1)-El)-R*(Isyn(i)))*(delT/tau);  
            v(i)=v(i-1)+delV;    
            if (spkOn==1)
                if (v(i)>thetaval)
                    v(i)=0;
                    v(i+1)=El;
                    i=i+1;
                end
            end
            i=i+1;
        end 
    end
