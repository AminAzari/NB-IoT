%plot orig
clc
close all
clear all

% figure(1)
% plot(dRg,AUd ,'--b')
% hold on 
% plot(dRg,ADd ,'--r')
% figure(2) 
% plot(dRg,ALd/4 ,'--b')
% hold on

%  
% KK=107;
% for io=1:KK
%     [io,1]
%     name=['data',num2str(io)];
%     load(name)
%     Mm(io)=mLd;
% end
% MM=min(Mm);
% for io=1:KK
%     [io,2]
%     name=['data',num2str(io)];
%     load(name)
%     nnDd=Dd(:,1:MM,:);
%     NDd(:,:,Nit*(io-1)+1:Nit*io)=nnDd;
% end
% 
% for io=1:KK
%     [io,3]
%     name=['data',num2str(io)];
%     load(name)
%     Mm(io)=mL;
% end
% MMu=min(Mm);
% for io=1:KK
%     [io,4]
%     name=['data',num2str(io)];
%     load(name)
%     nnDu=Du(:,1:MMu,:);
%     NDu(:,:,Nit*(io-1)+1:Nit*io)=nnDu;
% end
% 
% 
% for io=1:KK
%     [io,5]
%     name=['data',num2str(io)];
%     load(name)
%     Mm(io)=mLc;
% end
% MMc=min(Mm);
% for io=1:KK
%     [io,6]
%     name=['data',num2str(io)];
%     load(name)
%     nnDc=Dc(:,1:MMc,:);
%     NDc(:,:,Nit*(io-1)+1:Nit*io)=nnDc;
% end
% DDc=mean(mean(NDc,3),2);
% 
% % for io=1:KK
% %     [io,7]
% %     name=['data',num2str(io)];
% %     load(name)
% %     Mm(io)=mLtu;
% % end
% % MMtu=min(Mm);
% % for io=1:KK
% %     [io,8]
% %     name=['data',num2str(io)];
% %     load(name)
% %     nnDtu=Dtu(:,1:MMtu,:);
% %     NDtu(:,:,Nit*(io-1)+1:Nit*io)=nnDtu;
% % end
% % DDtu=mean(mean(NDtu,3),2);
% 
% % for io=1:KK
% %     [io,9]
% %     name=['data',num2str(io)];
% %     load(name)
% %     Mm(io)=mLrd;
% % end
% % MMrd=min(Mm);
% % for io=1:KK 
% %     [io,10]
% %     name=['data',num2str(io)];
% %     load(name)
% %     nnDrd=Drd(:,1:MMrd,:);
% %     NDrd(:,:,Nit*(io-1)+1:Nit*io)=nnDrd;
% % end
% % DDrd=mean(mean(NDrd,3),2);
% % % 
%  save data_KK107.mat

 %% run
 load datatestm.mat
ADt=sDd1(8,:);
AUt=sDu1(8,:);
ALt=sLT1(8,:);
% ADd=sDd1(:,16);
% AUd=sDu1(:,16);
% ALd=sLT1(:,16);

% figure(1)
% plot(tRg,AUt ,'--b')
% hold on 
% plot(tRg,ADt ,'--r')
% figure(2) 
% plot(tRg,ALt ,'--b')
% hold on

load data_KK107.mat
bDt=(f1*c1+f2*c2)*u; 
mQ=(Gu+Gd)*tRg+lab*d;
DraR=0.5*d+0.5*mQ*bDt+c*u;
DDu=mean(mean(NDu,3),2);  
DDd=mean(mean(NDd,3),2); 

Dsy=Dsy1;
Esyn=Pl*Dsy; 
Erar=Pl*DraR;%(d/2+1*0.002);
Era=PI*tRg/2+c*tau*(Pc+al*Pt);
Etx=PI*(DDu-c*ut/1000)+c*(Pc+al*Pt)*ut/1000;
Erx=PI*(DDd-c*dt/1000)+c*Pl*dt/1000;
Eu=Esyn+(Era+Erar)+Etx';
Ed=Esyn+(Era+Erar)+Erx';
E0=1000;  
   
LT=E0./(12*p*Eu+12*(1-p)*Ed);
 
figure(2) 
yyaxis right
plot(tRga,AUt ,'--b')
hold on  
plot(tRg,DDu+Dsy1+c*tau,'b') 
plot(tRga,ADt ,'--r')
plot(tRg,DDd+Dsy1+c*tau,'r')
hold off
yyaxis left
hold on 
plot(tRga,ALt,'--g')
hold on
plot(tRg,LT ,'g')
hold off
legend('An-U','U','An-D','D','A-L','L')


% figure(2)
% plot(tRg(1:1:end),LT(1:1:end),'r')
% hold off
% legend('An-L','Si-L')

