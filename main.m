% IHN
% clc
% clear all
% close all
function main(Nit,ind)
% Nit=1;
% ind=1000;
param
% Nit=5;
co=0;
% tRg=[0.002, 0.003,0.004,0.005,0.006,0.007, 0.008,0.01:0.005:   .2 ];
% tRg=[0.025:0.005:.2 ];
tRg=0.2;
for t=tRg
    co=co+1;
    for it=1:Nit
        %         t=0.028;
        d=0.01;
        % 0.5,1.5,2.5  NPRACH in ul
        %-0.5 ref in dl;-1.5 bs initiated control
        % 1:N--> device in ul; % j1:jN--> cont device in ul;
        %-1:-N--> device in dl; %-j1:-jN--> cont device in dl;
        T_max=1000*1000;
        
        GAc_t=zeros(1,T_max);
        GAc_t(t*1000:t*1000:T_max)=ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul=zeros(1,T_max);
        Ac_ul(t*1000:t*1000:T_max)=0.5*ones(1,length([t*1000:t*1000:T_max]));
        
        Ac_ul(t*1000+1:t*1000:T_max+1)=1.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+2:t*1000:T_max+2)=2.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+3:t*1000:T_max+3)=3.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+4:t*1000:T_max+4)=4.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+5:t*1000:T_max+5)=5.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+6:t*1000:T_max+6)=6.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+7:t*1000:T_max+7)=7.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+8:t*1000:T_max+8)=8.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+9:t*1000:T_max+9)=9.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+10:t*1000:T_max+10)=10.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+11:t*1000:T_max+11)=11.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+12:t*1000:T_max+12)=12.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+13:t*1000:T_max+13)=13.5*ones(1,length([t*1000:t*1000:T_max]));
        Ac_ul(t*1000+14:t*1000:T_max+14)=14.5*ones(1,length([t*1000:t*1000:T_max]));
        
        GAc_d=zeros(1,T_max); 
        GAc_d(1+d*1000:d*1000:T_max)=ones(1,length([d*1000+1:d*1000:T_max]));
        Ac_dl=zeros(1,T_max);
        Ac_dl(1:bt*1000:T_max)=-0.5*ones(1,length([1:bt*1000:T_max]));
        Ac_dl(2:bt*1000:T_max)=-0.5*ones(1,length([2:bt*1000:T_max]));
        
        for T=2:T_max/2
            if (GAc_d(T)==1)
                Kb=poissrnd(lambdab*d,1,1);
                if(Kb>0)
                    for jc=1:Kb
                        for ic=1:1000
                            if( sum(Ac_dl(T+ic:T+ic+u*c*1000-1))==0)
                                fT=T+ic;
                                break
                            end
                        end
                        Ac_dl(fT:fT+u*c*1000-1)=-1.5*ones(1,u*c*1000);
                    end
                end
            end
        end
        jjd=1;
        jju=1;
        jjc=1;
        jjd2=1;
        jju2=1;
        jjc2=1;
        
        er=0;
        erc=0;
        er2=0;
        erc2=0;
        
        Et1=zeros(1,T_max/2+1);
        Et2=zeros(1,T_max/2+1);
        Ev1=zeros(3,floor(N*S/(24*3600)*T_max/1000));
        Ev2=zeros(3,floor(N*S/(24*3600)*T_max/1000));
        DuEv1=zeros(10,floor(N*S/(24*3600)*T_max/1000));
        DuEv2=zeros(10,floor(N*S/(24*3600)*T_max/1000));
        DdEv1=zeros(10,floor(N*S/(24*3600)*T_max/1000));
        DdEv2=zeros(10,floor(N*S/(24*3600)*T_max/1000));
        lTim=zeros(N,floor(N*S/(24*3600)*T_max/1000));
        lTim2=zeros(N,floor(N*S/(24*3600)*T_max/1000));
        e1=1;e2=1;
        
                
        for T=2:T_max/2
            if(mod(T,T_max/2/10)==0)
                [co, it,T/(T_max/2),1,erc,erc2 ]
            end
            Ku=0;Kd=0;Kun=0;Kdn=0;
            Ku2=0;Kd2=0;Kun2=0;Kdn2=0;
            
            if (GAc_t(T-1)==1)
                Ku=poissrnd(N*S/(24*3600)*p*t*f1,1,1);
                Kd=poissrnd(N*S/(24*3600)*(1-p)*t*f1,1,1);
                Ku2=poissrnd(N*S/(24*3600)*p*t*f2,1,1);
                Kd2=poissrnd(N*S/(24*3600)*(1-p)*t*f2,1,1);
            end
            if(Ku+Kd>0)
                Sn=zeros(1,Ku+Kd);
                for ii=1:Ku+Kd
                    Sn(ii)=randperm(M,1);
                end
                Sc=unique(Sn);
                ps=length(Sc)/(Ku+Kd);
                Ku=ceil(Ku*ps);
                Kd=ceil(Kd*ps);
                Kun=randperm(N,Ku);
                Kdn=randperm(N,Kd);
                
            end
            if(Ku2+Kd2>0)
                Sn=zeros(1,Ku2+Kd2);
                for ii=1:Ku2+Kd2
                    Sn(ii)=randperm(M,1);
                end
                Sc=unique(Sn);
                ps=length(Sc)/(Ku2+Kd2);
                Ku2=ceil(Ku2*ps);
                Kd2=ceil(Kd2*ps);
                Kun2=randperm(N,Ku2);
                Kdn2=randperm(N,Kd2);
            end
            if(Ku+Kd>0)
                Et1(T)=1;
                Ev1(:,e1)=[T;Ku;Kd];
                DuEv1(1:Ku,e1)=Kun';
                DdEv1(1:Kd,e1)=Kdn';
                e1=e1+1;
                
                for jc=1:Ku
                    fT=0;
                    for ic=1:100000
                        if( sum(Ac_dl(T+ic:T+ic+u*c*1000-1))==0)
                            if(GAc_d(T+ic)==1 || Ac_dl(T+ic-1)==-1.5 || ic>d)
                                fT=T+ic;
                                break
                            end
                        end
                    end 
                    if(fT>0)
                        Ac_dl(fT:fT+ u*c*1000-1)=1j*Kun(jc)*ones(1, u*c*1000);
                        Dc(co,jjc,it)=(fT+ u*c*1000-1-T)/1000;%+t/2;
                        jjc=jjc+1;
                        lTim(Kun(jc),e1-1)=fT+ u*c*1000-1+1;
                    else
                        erc=erc+1;
                    end
                end
                
                for jc=1:Kd
                    fT=0;
                    for ic=1:100000
                        if( sum(Ac_dl(T+ic:T+ic+u*c*1000-1))==0)
                            if(GAc_d(T+ic)==1 || Ac_dl(T+ic-1)==-1.5  )
                                fT=T+ic;
                                break
                            end
                        end
                    end
                    if(fT>0)
                        Ac_dl(fT:fT+ u*c*1000-1)=-1j*Kdn(jc)*ones(1, u*c*1000);
                        Dc(co,jjc,it)=(fT+ u*c*1000-1-T)/1000;%+t/2;
                        jjc=jjc+1;
                        
                    else
                        erc =erc+1;
                    end
                    
                end
                
                
            end
            
            
            if(Ku2+Kd2>0)
                Et2(T)=1;
                Ev2(:,e2)=[T;Ku2;Kd2];
                DuEv2(1:Ku2,e2)=Kun2';
                DdEv2(1:Kd2,e2)=Kdn2';
                e2=e2+1;
                for jc=1:Ku2
                    fT=0;
                    for ic=1:100000
                        if( sum(Ac_dl(T+ic:T+ic+u*c2*1000-1))==0)
                            if(GAc_d(T+ic)==1 || Ac_dl(T+ic-1)==-1.5|| ic>2*d)
                                fT=T+ic;
                                break
                            end
                        end
                    end
                    if(fT>0)
                        Ac_dl(fT:fT+ u*c2*1000-1)=1j*Kun2(jc)*ones(1, u*c2*1000);
                        Dc2(co,jjc,it)=(fT+ u*c2*1000-1-T)/1000;%+t/2;
                        jjc2=jjc2+1;
                        lTim2(Kun2(jc),e2-1)=fT+ u*c2*1000-1+1;
                    else
                        erc2 =erc2 +1;
                    end
                end
                
                for jc=1:Kd2
                    fT=0;
                    for ic=1:100000
                        if( sum(Ac_dl(T+ic:T+ic+u*c2*1000-1))==0)
                            if(GAc_d(T+ic)==1 || Ac_dl(T+ic-1)==-1.5 || ic>2*d )
                                fT=T+ic;
                                break
                            end
                        end
                    end
                    if(fT>0)
                        Ac_dl(fT:fT+ u*c2*1000-1)=-1j*Kdn2(jc)*ones(1, u*c2*1000);
                        Dc2(co,jjc2,it)=(fT+ u*c2*1000-1-T)/1000;%+t/2;
                        jjc2=jjc2+1;
                        
                    else
                        erc2 =erc2 +1;
                    end
                    
                end
                
            end
        end
        
        
        e1=0;e2=0;
        
        for T=2:T_max/2
            if(mod(T,T_max/2/10)==0)
                [co, it,T/(T_max/2),2 ,er,er2 ]
            end
            Ku=0;Kd=0;Kun=0;Kdn=0;
            Ku2=0;Kd2=0;Kun2=0;Kdn2=0;
            if (Et1(T)==1)
                e1=e1+1;
                Ku=Ev1(2,e1);
                Kd=Ev1(3,e1);
                Kun=DuEv1(1:Ku,e1);
                Kdn=DdEv1(1:Kd,e1);
            end
            
            if (Et2(T)==1)
                e2=e2+1;
                Ku2=Ev2(2,e2);
                Kd2=Ev2(3,e2);
                Kun2=DuEv2(1:Ku2,e2);
                Kdn2=DdEv2(1:Kd2,e2);
            end
            
            
            if(Ku+Kd>0)
                
                fTm=0;
                for jc=1:Ku
                    fTo=max(T,lTim(Kun(jc),e1));
                    fT=0;
                    for ijc=1:ut
                        
                        if(ijc>1 && fTm==0 && fT>0)
                            fTm=fT;
                        end
                        if(ijc>1)
                            fTo=fT;
                        end
                        for ic=1:100000
                            if( sum(Ac_ul(fTo+ic:fTo+ic+TTI*c*1000-1))==0)
                                fT=fTo+ic;
                                break
                            end
                        end
                        if(fT>0)
                            Ac_ul(fT:fT+ TTI*c*1000-1)=Kun(jc)*ones(1, TTI*c*1000);
                            
                        else
                            er =er+1;
                            break
                        end
                    end
                    if(fT>0) 
                        Du(co,jju,it)=(fT+ TTI*c*1000-1-T)/1000+t/2;
                        Dtu(co,jju,it)=(fT+ TTI*c*1000-1-fTm)/1000;
                        jju=jju+1; 
                    end
                end
                
                
                
                for jc=1:Kd
                    fTo=T;
                    fT=0;
                    for ijc=1:dt
                        if(ijc>1 && fTm==0 && fT>0)
                            fTm=fT;
                        end
                        if(ijc>1)
                            fTo=fT;
                        end
                        for ic=1:100000
                            if( sum(Ac_dl(fTo+ic:fTo+ic+TTI*c*1000-1))==0)
                                fT=fTo+ic;
                                break 
                            end
                        end
                        if(fT>0)
                            Ac_dl(fT:fT+ TTI*c*1000-1)=-Kdn(jc)*ones(1, TTI*c*1000);
                        else
                            er=er+1;
                            break 
                        end
                    end  
                    if(fT>0) 
                        Dd(co,jjd,it)=(fT+ TTI*c*1000-1-T)/1000+t/2; 
                        Drd(co,jjd,it)=(fT+ TTI*c*1000-1-fTm)/1000; 
                        jjd=jjd+1;   
                    end
                    
                end
            end
            
            
            if(Ku2+Kd2>0)
                fTm=0;
                for jc=1:Ku2
                    fTo=max(T,lTim2(Kun2(jc),e2));
                    fT=0;
                    for ijc=1:ut
                        
                        if(ijc>1 && fTm==0 && fT>0)
                            fTm=fT;
                        end
                        if(ijc>1)
                            fTo=fT;
                        end
                        for ic=1:100000
                            if( sum(Ac_ul(fTo+ic:fTo+ic+TTI*c2*1000-1))==0)
                                fT=fTo+ic;
                                break
                            end
                        end
                        if(fT>0)
                            Ac_ul(fT:fT+ TTI*c2*1000-1)=Kun2(jc)*ones(1, TTI*c2*1000);
                        else
                            er2=er2+1;
                            break
                        end
                    end
                    if(fT>0)
                        Du2(co,jju2,it)=(fT+ TTI*c2*1000-1-T)/1000+t/2;
                        Dtu2(co,jju2,it)=(fT+ TTI*c2*1000-1-fTm)/1000;
                        jju2=jju2+1;
                    end
                end
                
                
                
                for jc=1:Kd2
                    fTo=T;
                    fT=0;
                    for ijc=1:dt
                        if(ijc>1 && fTm==0 && fT>0)
                            fTm=fT;
                        end
                        if(ijc>1)
                            fTo=fT;
                        end
                        for ic=1:100000
                            if( sum(Ac_dl(fTo+ic:fTo+ic+TTI*c2*1000-1))==0)
                                fT=fTo+ic;
                                break
                            end
                        end
                        if(fT>0)
                            Ac_dl(fT:fT+ TTI*c2*1000-1)=-Kdn2(jc)*ones(1, TTI*c2*1000);
                        else
                            er2=er2+1;
                            break
                        end
                    end
                    if(fT>0)
                        Dd2(co,jjd2,it)=(fT+ TTI*c2*1000-1-T)/1000+t/2;
                        Drd2(co,jjd2,it)=(fT+ TTI*c2*1000-1-fTm)/1000;
                        jjd2=jjd2+1;
                    end
                end
            end
            
            
        end
        
         
    end
end
% tRg=dRg;
[m1,m2,m3]=size(Du);
Lgt=zeros(1,Nit);
for j=1:Nit
    for l=1:m2
        for i=1:length(tRg)
            if(Du(i,l,j)==0)
                Lgt(j)=l-1;
                break
            end
        end
        if(Lgt(j)>0)
            break
        end
    end
    
end
mL=min(Lgt);
% nDu=Du(:,10:mL,:);

[m1,m2,m3]=size(Du2);
Lgt=zeros(1,Nit);
for j=1:Nit
    for l=1:m2
        for i=1:length(tRg)
            if(Du2(i,l,j)==0)
                Lgt(j)=l-1;
                break
            end
        end
        if(Lgt(j)>0)
            break
        end
    end
    
end
mL2=min(Lgt);
% nDu=Du(:,10:mL,:);



[m1,m2,m3]=size(Dd);
Lgt=zeros(1,Nit);
for j=1:Nit
    for l=1:m2
        for i=1:length(tRg)
            if(Dd(i,l,j)==0)
                Lgt(j)=l-1;
                break
            end
        end
        if(Lgt(j)>0)
            break
        end
    end
end
mLd=min(Lgt);
% nDd=Dd(:,3:mLd,:);


[m1,m2,m3]=size(Dd2);
Lgt=zeros(1,Nit);
for j=1:Nit
    for l=1:m2
        for i=1:length(tRg)
            if(Dd2(i,l,j)==0)
                Lgt(j)=l-1;
                break
            end
        end
        if(Lgt(j)>0)
            break
        end
    end
end
mLd2=min(Lgt);
% nDd=Dd(:,3:mLd,:);


[m1,m2,m3]=size(Dc);
Lgt=zeros(1,Nit);
for j=1:Nit
    for l=1:m2
        for i=1:length(tRg)
            if(Dc(i,l,j)==0)
                Lgt(j)=l-1;
                break
            end
        end
        if(Lgt(j)>0)
            break
        end
    end
end
mLc=min(Lgt);
% nDc=Dc(:,3:mLc,:);

[m1,m2,m3]=size(Dc2);
Lgt=zeros(1,Nit);
for j=1:Nit
    for l=1:m2
        for i=1:length(tRg)
            if(Dc2(i,l,j)==0)
                Lgt(j)=l-1;
                break
            end
        end
        if(Lgt(j)>0)
            break
        end
    end
end
mLc2=min(Lgt);
% nDc=Dc(:,3:mLc,:);



[m1,m2,m3]=size(Dtu);
Lgt=zeros(1,Nit);
for j=1:Nit
    for l=1:m2
        for i=1:length(tRg)
            if(Dtu(i,l,j)==0)
                Lgt(j)=l-1;
                break
            end
        end
        if(Lgt(j)>0)
            break
        end
    end
end
mLtu=min(Lgt);
% nDtu=Dtu(:,3:mLtu,:);


[m1,m2,m3]=size(Dtu2);
Lgt=zeros(1,Nit); 
for j=1:Nit
    for l=1:m2
        for i=1:length(tRg)
            if(Dtu2(i,l,j)==0)
                Lgt(j)=l-1;
                break
            end
        end
        if(Lgt(j)>0)
            break
        end
    end
end
mLtu2=min(Lgt);
% nDtu=Dtu(:,3:mLtu,:);


[m1,m2,m3]=size(Drd);
Lgt=zeros(1,Nit);
for j=1:Nit
    for l=1:m2
        for i=1:length(tRg)
            if(Drd(i,l,j)==0)
                Lgt(j)=l-1;
                break
            end
        end
        if(Lgt(j)>0)
            break
        end
    end
end
mLrd=min(Lgt);
% nDrd=Dtu(:,3:mLrd,:);

[m1,m2,m3]=size(Drd2);
Lgt=zeros(1,Nit);
for j=1:Nit
    for l=1:m2
        for i=1:length(tRg)
            if(Drd2(i,l,j)==0)
                Lgt(j)=l-1;
                break
            end
        end
        if(Lgt(j)>0)
            break
        end
    end
end
mLrd2=min(Lgt);
% nDrd=Dtu(:,3:mLrd,:);


name=['data',num2str(ind)];
save(name)




% DDu=mean(mean(nDu,3),2);
% DDd=mean(mean(nDd,3),2);
% plot(tRg,DDu,'b')
% hold on
% plot(tRg,DDd,'r')
% hold off

% for i=1:10
%     plot(mean(nDd(i,:,:),3))
%     hold on
% end

