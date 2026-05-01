% S_DEMO - demonstration of capabilities.
% Short example of capabilities of S_TIDE.
%%%%== https://www.researchgate.net/profile/Haidong-Pan ================!!!!
%%%%========== https://space.bilibili.com/548420599     ================!!!
%%%%========== http://blog.sina.com.cn/panhaidongphd    ================!!!
%%
% demonstration of s_calculate_coefficient.m
x=[1 2 3 4 5];y=[0 1 0 -1 0];
xi=1:0.05:5;
IDP1=5;
nobs=length(xi);
S1 = s_calculate_coefficient(IDP1,nobs); %use spline interpolation 
S2 = l_calculate_coefficient(IDP1,nobs);  %use linear interpolation
S3 = sinc_calculate_coefficient(IDP1,nobs);  %use sinc interpolation
yi_spline=S1*y';
yi_linear=S2*y';
yi_sinc=S3*y';
plot(x,y,'ro','MarkerSize',8.0);hold on;plot(xi,yi_spline,'b','linewidth',1.1);
plot(xi,yi_linear,'m','linewidth',1.1);plot(xi,yi_sinc,'k','linewidth',1.1)
set(gca,'Fontsize',15,'linewidth',1.2);
legend('Given','Spline','Linear','Sinc')
ylim([-1.2 1.2])

%%
% demonstration of s_tide.m
load kushiro.mat   % hourly elevation data at Kushiro (Japan) from 1993-01-01 to 2012-12-31, obtained from PSMSL
load imf.mat      % the modes obtained by Empirical Mode Decomposition (EMD) of the kushiro elevation data.
% If IDP is equal to 2, s_tide will calculate the linear trend of water level!
%If IDP is equal to 1, s_tide will calcualte a constant!
figure();
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro,1,5,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'Ssa';'Sa'},10,1,'spline','ols');
plot(St);hold on
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro,2,5,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'Ssa';'Sa'},10,1,'spline','ols');
plot(St,'r')
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro,3,5,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'Ssa';'Sa'},10,1,'spline','ols');
plot(St,'k')
legend('IP=1','IP=2','IP=3')
%calculate given frequencies(non-tidal)
figure();
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro,2,21,[1/36,1/72],2,1,'spline','ols');%period(36 hours and 72 hours)
plot(Ht(1,:),'r'); hold on;plot(Ht(2,:),'k')

%
figure();
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro,21,21,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'Ssa';'Sa'},10,1,'spline','ols');
plot(St,'k','linewidth',1.2)
hold on
plot(imf(16,:)+imf(15,:)+imf(14,:)+imf(13,:)+imf(12,:),'r','linewidth',1.2)
legend('S\_TIDE','EMD')
%legend('boxoff')
set(gca,'Fontsize',15,'linewidth',1.2)
set(gca,'XTick',0:8760*2:8760*20)
set(gca,'XTickLabel',{'1993','1995','1997','1999','2001','2003','2005','2007','2009','2011','2013'})
ylabel('Sea Level/mm','Fontsize',15)
xlabel('Year','Fontsize',15)
grid on
xlim([0 8760*20])
%plot the amplitude and phase and their 95% confidence intervals
for i=1:2
figure()
plot(Ht(i,:));hold on
plot(Ht(i,:)+Htint(i,:),'r')
plot(Ht(i,:)-Htint(i,:),'r')
end
for i=1:2
figure()
plot(Gt(i,:));hold on
plot(Gt(i,:)+Gtint(i,:),'r')  
plot(Gt(i,:)-Gtint(i,:),'r')
end
%%
% demonstration of s_nodal_correction.m
% You must have installed T_TIDE to run the code below
[NAME,FREQ,TIDECON,XOUT]=t_tide(kushiro(1:8760*3),'interval',1,'rayleigh',['M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2']); %without nodal correction
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro(1:8760*3),4,4,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'},8,1,'spline','ols');
subplot(2,1,1);plot(Ht(1,:),'r');hold on;plot(TIDECON(6,1)*ones(1,8760*3),'k');ylabel('Amplitude'); title('nodal uncorrected')
subplot(2,1,2);plot(Gt(1,:),'r');hold on;plot(TIDECON(6,3)*ones(1,8760*3),'k');ylabel('Phase')
legend('S\_TIDE','T\_TIDE')

[NAME,FREQ,TIDECON,XOUT]=t_tide(kushiro(1:8760*3),'interval',1,'rayleigh',['M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'],'latitude',42.5,'start time',[1993,01,01,00]);
[Hc,Gc,ff,uu,vv]=s_nodal_correction(Ht,Gt,[1993,01,01],8760*3,1,42.5,ju); figure()
subplot(2,1,1);plot(Hc(1,:),'r');hold on;plot(TIDECON(6,1)*ones(1,8760*3),'k');ylabel('M2 Amplitude');title('nodal corrected')
subplot(2,1,2);plot(Gc(1,:),'r');hold on;plot(TIDECON(6,3)*ones(1,8760*3),'k');ylabel('M2 Phase')
legend('S\_TIDE','T\_TIDE')

figure()
subplot(2,1,1);plot(Hc(3,:),'r');hold on;plot(TIDECON(4,1)*ones(1,8760*3),'k');ylabel('K1 Amplitude');title('nodal corrected')
subplot(2,1,2);plot(Gc(3,:),'r');hold on;plot(TIDECON(4,3)*ones(1,8760*3),'k');ylabel('K1 Phase')
legend('S\_TIDE','T\_TIDE')

figure()
subplot(2,1,1);plot(Hc(7,:),'r');hold on;plot(TIDECON(3,1)*ones(1,8760*3),'k');ylabel('P1 Amplitude');title('nodal corrected')
subplot(2,1,2);plot(Gc(7,:),'r');hold on;plot(TIDECON(3,3)*ones(1,8760*3),'k');ylabel('P1 Phase')
legend('S\_TIDE','T\_TIDE')
%%
%extract 18.61 year cycle from monthly K1 tidal amplitudes using S_TIDE
load K1_amp_hilo.mat    
% the IP number for 18.61 year cycle is one which means the amplitude and phase of 18.61 year cycle are constant
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(K1_amp,2,1,[1/(18.61*8760)],1,720,'spline','ols');
figure();plot(K1_amp,'r','linewidth',1.2);hold on;plot(xout,'k','linewidth',1.2)
 mean(Ht)/mean(St)*100  % theoretical value is 11.5%, but the actual value is 11.72%

% If the IP number for 18.61 year cycle is zero, s_tide will only calculate the linear trend
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(K1_amp,2,0,[1/(18.61*8760)],1,720,'spline','ols');
figure();plot(K1_amp,'r','linewidth',1.2);hold on;plot(xout,'k','linewidth',1.2)
%%
% OLS vs robustfit, this example was first added in S_TIDE v1.14
load kushiro.mat   % hourly elevation data at Kushiro (Japan) from 1993-01-01 to 2012-12-31, obtained from PSMSL
load imf.mat 
tic
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro(1:8760*9),10,10,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'},8,1,'spline','ols');
toc
%ols needs 0.874286s
tic
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro(1:8760*9),10,10,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'},8,1,'spline','robustfit');
toc
%robustfit needs 20.675073s
% Ols is much faster than robustfit, while their results are nearly same!

%%
%demonstration of s_tide_m1.m
%this example was first added in S_TIDE v1.16
%Reference:Wang, D., H. Pan, G. Jin and X. Lv (2020), Seasonal variation of the main tidal constituents in the Bohai Bay, Ocean Science.
load kushiro.mat 
[St,Ht,Gt,H,G,coef,xout,ju,coefint]=s_tide_m1(kushiro,200,50,{'M2';'K1';'O1';'S2';'N2';'K2';'P1';'Q1'},8,2,1);
plot(Ht(1,:),'k');hold on;plot(Ht(2,:),'r')
legend('M2','K1');xlabel('Time(hour)');ylabel('Amplitude(mm)')
%ÎÒÃÇ¿ÉÒÔÇå³þµÄ¿´µ½M2ºÍK1·Ö³±µÄÄê±ä»¯ºÍ18.61Äê±ä»¯
%We can clearly see the annual and nodal variations in M2 and K1 amplitudes

%%
%demonstration of s_tide_m2.m (the function s_tide_m1.m can be seen as a special case of s_tide_m2.m)
%this example was first added in S_TIDE v1.16
load kushiro.mat 
[St,Ht,Gt,coef,xout,ju,coefint,Stint,Htint,Gtint,consti]=s_tide_m2(kushiro,20,50,5,{'M2';'K1';'O1';'S2';'N2';'K2';'P1';'Q1'},8,2,1);
figure();plot(Ht(1,:),'k');hold on;plot(Ht(2,:),'r');plot(Ht(3,:),'b')
legend('M2','K1','O1');xlabel('Time(hour)');ylabel('Amplitude(mm)')
%ÎÒÃÇ¿ÉÒÔÇå³þµÄ¿´µ½M2ºÍK1·Ö³±µÄÄê±ä»¯ºÍ18.61Äê±ä»¯
%We can clearly see the annual and nodal variations in M2 and K1 amplitudes

%%
%demonstration of s_tide_m3.m
%this example was first added in S_TIDE v1.16
load('satellite.mat');
 cons={'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'SA';'SSA'};lat=14.1633;
 stime=[1992,10,31,19,43,37.71];
 [St1,Ht1,Gt1,coef,xout1,ju,Stint1,Htint1,Gtint1,aa1,bb1,sigma_aa1,sigma_bb1]=s_tide_m3(satellite,2,1,cons,10,9.9156*24,'spline','ols',lat,stime);
  
 plot(satellite,'*-','linewidth',1.1);hold on;plot(xout1,'r','linewidth',1.1);
 set(gca,'Fontsize',12,'linewidth',1.2)
 xlim([1 101]);set(gca,'XTick',[1 51 101])
 ylabel('Sea Level(m)');xlabel('Time')
  set(gca,'XTickLabel',{'1992/10/31','1994/03/10', '1995/07/19'})
    legend('Observation','Hindcast')
%%
% this example was first added in S_TIDE v1.17
load kushiro.mat %data start from 1993/01/01
[nameu,fu,tidecon,xout]=t_tide(kushiro(1:8767)); % without nodal corrections
[nameu,fu,tidecon,xout]=t_tide(kushiro(1:8767),'latitude',42.5,'start time',[1993,01,01,00]);%with nodal corrections

% S_TIDE results without nodal corrections
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint]=s_tide(kushiro(1:8767),1,1,'autoselected','autoselected',1,'spline','ols');
%compared to T_TIDE,the results are nearly same
nobs=length(kushiro(1:8767));nobsu=nobs-rem(nobs-1,2);% makes series odd to give a center point
ctime=datenum(1993,1,1)+(ceil(nobsu/2)-1)/24.0; % using center time for nodal correction
[v,u,f]=t_vuf('nodal',ctime,ju,42.5);
u=360*u;%convert to degree
v=360*v;
Hc=Ht(:,1)./f; %nodal corrected amplitudes
Gc=Gt(:,1)+u+v; %nodal corrected phases
%%
%demonstration of s_tide_m4.m
%this example was first added in S_TIDE v1.17_update
load('satellite.mat');
 cons={'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'};lat=14.1633;
 stime=[1992,10,31,19,43,37.71];
[St,Ht,Gt,coef,xout,ju,coefint,Stint,Htint,Gtint]=s_tide_m4(satellite,15,6,1,cons,8,4,9.9156*24,lat,stime);

 plot(satellite,'*-','linewidth',1.1);hold on;plot(xout,'r','linewidth',1.1);
 set(gca,'Fontsize',12,'linewidth',1.2)
 xlim([1 101]);set(gca,'XTick',[1 51 101])
 ylabel('Sea Level(m)');xlabel('Time')
  set(gca,'XTickLabel',{'1992/10/31','1994/03/10', '1995/07/19'})
    legend('Observation','Hindcast')
%%
%demonstration of s_tide_m5.m
%this example was first added in S_TIDE v1.18
sl=ncread('h572a.nc','sea_level');time1=ncread('h572a.nc','time');
time=datevec(datenum([1800 01 01])+time1);
time(1,:)  %observations at Astoria start from 1925/01/25

[St,Ht,Gt,coef,xout,ju,coefint,Stint,Htint,Gtint,consti]=s_tide_m5(sl(1:8767*20),20,15,5,{'M2';'K1';'O1';'S2';'N2';'K2';'P1';'Q1'},8,2,time(1:8767*20,:));

%%
%demonstration of s_tide_m6.m
%this example was first added in S_TIDE v1.18
%¼ÆËãµ÷ºÍ³£Êý calculate harmonic constants
clear
load('satellite.mat');
 cons={'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'};lat=14.1633;
 stime=[1992,10,31,19,43,37.71];
 for i=1:length(satellite)
     time(i,:)=datevec(datenum(stime)+(i-1)*9.9156);
 end
 satellite(2)=[];time(2,:)=[];
 
[St2,Ht2,Gt2,coef,xout,ju,coefint,Stint,Htint,Gtint,consti]=s_tide_m6(satellite,1,1,1,cons,8,1,lat,time);

%%
%demonstration of s_tide for tidal currents
%this example was first added in S_TIDE v 1.19
%³±Á÷ÍÖÔ²»æÖÆºÍÍÖÔ²²ÎÊý¼ÆËãÊµÀý
load tidalcurrents.mat

[Stv,Htv,Gtv,coef,xoutv,ju,Stint,Htint,Gtint,aa,bb,namev,fv,constiv]=s_tide_m8(v,1,1,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'},8,1,'spline','robustfit',20.5,[2010,08,23,01,00,00],'corrected');
[Stu,Htu,Gtu,coef,xoutu,ju,Stint,Htint,Gtint,aa,bb,nameu,fu,constiu]=s_tide_m8(u,1,1,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'},8,1,'spline','robustfit',20.5,[2010,08,23,01,00,00],'corrected');
% 20.5 is the latitude of the observation,[2010,08,23,01,00,00] is the start time of observation, users only need to change these two inputs when perform their own harmonic analysis.
% v and u are observed tidal currents at y and x axis
s_plot_tidal_ellipse(Htu(1,1),Gtu(1,1),Htv(1,1),Gtv(1,1),'M2') %plot M2 tidal ellipse
s_plot_tidal_ellipse(Htu(3,1),Gtu(3,1),Htv(3,1),Gtv(3,1),'K1') %plot K1 tidal ellipse
s_plot_tidal_ellipse(Htu(5,1),Gtu(5,1),Htv(5,1),Gtv(5,1),'N2') %plot N2 tidal ellipse
%calculate the parameters of M2 tidal ellipse
[maxc,minc,maxt,mint,maxp,minp,rr]=s_estimate_tidal_ellipse(Htu(1,1),Gtu(1,1),Htv(1,1),Gtv(1,1),'M2') 
% output:
%maxc: maximum tidal current speed (longer axis of tidal ellipse) ³±Á÷ÍÖÔ²µÄ³¤Öá
%minc: minimum tidal current speed (shorter axis of tidal ellipse) ³±Á÷ÍÖÔ²µÄ¶ÌÖá
%maxt: the time of maximum tidal current (start from 0 hour)
%mint: the time of minimum tidal current (start from 0 hour)
%maxp: the angle of maximum tidal current
%minp: the angle of minimum tidal current
%rr:   rotation rate(Ðý×ªÂÊ)which is the ratio of minc to maxc. 
%(The tidal flow vector rotates clockwise with time when rotation rate is negative) 
% this is a dynamic example which can help you understand the meaning of the parameters of tidal ellipse


%Ò»¸ö¶¯Ì¬ÊµÀý(ÔËÐÐ³ÌÐò»á¿´µ½³±Á÷ÍÖÔ²µÄÉú³É¹ý³Ì£©À´½âÊÍ³±Á÷ÍÖÔ²ÒªËØµÄ¾ßÌåº¬Òå
uam=38.6;uph=247; %sind(uph-vph)>0³±Á÷ÍÖÔ²¾ÍÊÇË³Ê±ÕëÐý×ª£¬sind(uph-vph)<0¾ÍÊÇÄæÊ±Õë
vam=25.0;vph=222;
tides={'M2'};
load t_constituents
ju=strmatch(upper(tides),const.name);
fu=const.freq(ju);
t=[0:1/12:36]; %from 0 hour to 36 hours ´Ó0Ê±¿Ìµ½36Ð¡Ê±
u=uam*cos(fu*t*2*pi-uph/180*pi);
v=vam*cos(fu*t*2*pi-vph/180*pi);
L=sqrt(u.^2+v.^2);
for i=1:160
subplot(2,1,1);plot(u(i),v(i),'k.');
title([ 'time=',num2str(t(i)),'hours'])
hold on;xlim([-40 40]);ylim([-30 30]);
xlabel('u(cm/s)','fontsize',12);ylabel('v(cm/s)','fontsize',12);
pause(0.04)
subplot(2,1,2);plot(t(i),L(i),'k.');hold on
xlim([t(1) t(160)]);ylim([-10 60])
xlabel('time(hours)','fontsize',12);ylabel('the speed of tidal current vector(cm/s)','fontsize',12)
pause(0.04)
end
[a1,a2]=findpeaks(L);maxt=t(a2);
[a3,a4]=findpeaks(-1*L);mint=t(a4);
maxc=max(L);minc=min(L);
line([t(a2(1)) t(a2(1))],[-10 maxc],'color','r');
%line([0 t(a2(1))],[maxc maxc],'color','r');
text(0.5,20,'maxt = 2.083 hour');
text(0.5,25,'maxc = 45.09')

line([t(a2(2)) t(a2(2))],[-10 maxc],'color','r');
%line([0 t(a2(2))],[maxc maxc],'color','r');
text(8.4,22,'maxt = 8.29 hour');
text(8.4,27,'maxc = 45.09')

line([t(a4(1)) t(a4(1))],[-10 minc],'color','r');
%line([0 t(a4(1))],[minc minc],'color','r');
text(4.5,22,'mint = 5.17 hour');
text(4.5,27,'minc = 9.05')

line([t(a4(2)) t(a4(2))],[-10 minc],'color','r');
%line([0 t(a4(2))],[minc minc],'color','r');
text(11,21,'mint = 11.38 hour');
text(11,26,'minc = 9.05')

subplot(2,1,1);line([-40 40],[0 0],'color','r');
line([0 0],[-30 30],'color','r');
line([u(a2(1)) u(a2(2))],[v(a2(1)) v(a2(2))]);
text(10,3,'31.76¡ã')
text(5,-13,'There are two angles for maximum tidal currents','Fontsize',12)
text(7,-18,'One is 31.76¡ã, the other is 31.76¡ã+ 180¡ã','Fontsize',12)

%%
%demonstration of s_consturct which can construct tidal water levels by given the information of tides
%this example was first added in S_TIDE v 1.19 update
%¸ù¾Ý¸ø¶¨µÄ·Ö³±Õñ·ùºÍ³Ù½Ç¹¹ÔìË®Î»
tides={'M2';'S2';'N2'};
amp=[1,0.5,0.2];pha=[180,0,90];dt=1;length=8760;
[sumconsti,consti]=s_construct(tides,amp,pha,dt,length);
plot(consti(1,:),'r.-')
xlim([1 48]);hold on
plot(consti(2,:),'k.-');plot(consti(3,:),'m.-')
legend('M2','S2','N2');xlabel('Time(hour)')
figure();plot(consti(1,:)+consti(2,:),'b-') %Ë®Î»°ëÔÂ±ä»¯
xlim([1 720*2]);xlabel('Time(hour)')
title('M2 tide + S2 tide') 
figure();plot(consti(1,:)+consti(3,:),'g-')%Ë®Î»ÔÂ±ä»¯
xlim([1 720*2]);xlabel('Time(hour)')
title('M2 tide + N2 tide')
figure();plot(consti(1,:).*consti(2,:),'b-') %Ë®Î»°ëÔÂ±ä»¯
xlim([1 720*2]);xlabel('Time(hour)')
title('M2 tide * S2 tide')
%plot S_TIDE toolbox logo  »æÖÆs_tide¹¤¾ß°üµÄlogo
figure();plot(sumconsti,'r-') 
xlim([1 720*2]);xlabel('Time(hour)')
text(300,0,'S\_TIDE','color','k','Fontsize',70)
set(gca,'Fontsize',13)
%%
%demonstration of s_nodal_cal which was added in S_TIDE v1.20
%¸ù¾Ý¸ø¶¨µÄ·Ö³±ºÍÊ±¼ä¼ÆËã½»µãÒò×ÓºÍ¶©Õý½Ç
[ff,uu,vv,N]=s_nodal_cal({'K1';'O1';'M2'},[1920,01,01],[2020,01,01],365,10);
subplot(2,1,1);
plot(ff(1,:),'k');hold on;plot(ff(2,:),'r');plot(ff(3,:),'b')
legend('K1','O1','M2');
set(gca,'XTick',0:25:2020)
set(gca,'XTickLabel',{'1920','1945','1970','1995','2020'})
title('(a) nodal factor')
subplot(2,1,2);
plot(uu(1,:),'k');hold on;plot(uu(2,:),'r');plot(uu(3,:),'b')
set(gca,'XTick',0:25:2020)
set(gca,'XTickLabel',{'1920','1945','1970','1995','2020'})
title('(b) nodal angle')
%%
clear;clc;
%Ê¹ÓÃkushiroÕ¾°Ë´óÖ÷Òª·Ö³±µ÷ºÍ³£ÊýÔ¤±¨Ë²Ê±³±Î»
% this example was first added in S_TIDE v1.20 update
tides={'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'};
dt=1;length=720;
amp=[288.23,133.23,250.80,193.39,35.1736,35.3202,82.4656,36.9610];
pha=[177.14,217.38,24.42,356.23,159.69,345.72,21.14,211.69];
[sumconsti2,consti2]=s_construct2(tides,amp,pha,dt,length,[1994,01,01],42.5);

load kushiro.mat  %the observed water levels at Kushiro start from 1993/01/01 00:00:00
plot(sumconsti2);hold on
plot(kushiro(8761:8761+720)-1807.48,'r')%1807.48 is mean water level
xlim([1 350])
%%
%demonstration of s_tide_m7.m
%this example was first added in S_TIDE v1.22   2021/08/16
load('satellite.mat');
 cons={'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2';'SA';'SSA'};lat=14.1633;
 stime=[1992,10,31,19,43,37.71];
 [St,Ht,Gt,coef,xout,ju,aa,bb,ta,tb,taint,tbint,coefint]=s_tide_m7(satellite,cons,10,9.9156*24,lat,stime);
  plot(aa(1,:));hold on;plot(bb(1,:),'k') %M2 tide 
%%
%demonstration of s_tide_m8.m
%this example was first added in S_TIDE v1.22  2021/08/16
load kushiro.mat
[St,Ht,Gt,coef,xout,ju,Stint,Htint,Gtint,aa,bb,nameu,fu,consti]=s_tide_m8(kushiro(1:8767),1,1,'autoselected','autoselected',1,'spline','robustfit',42.5,[1993,01,01,00,00,00],'corrected');

%%
%demonstration of s_tidalcharacter.m
%this example was first added in S_TIDE v1.22 update  2021/08/24
sl=ncread('h572a.nc','sea_level');time1=ncread('h572a.nc','time');
stime=datevec(datenum([1800 01 01])+time1(1,:))
[Htides,Htime,Ltide,Ltime,skewness]=s_tidalcharacter(sl(1:720),stime,1);

%%
%demonstration of s_alias.m
%this example was first added in S_TIDE v1.22 update2  2021/08/29
[Ta]=s_alias(9.9156,{'M2';'S2';'K1';'O1';'Msf'})
%T/P-Jason alias period: 62.11,58.74,173.19,45.71,30.19
[Ta]=s_alias(35,{'M2';'K1';'O1';'Mm'})
%Envisat alias period: 94.49,365.24,173.19,75.07,129.53

[Ta]=s_alias(9.9156,[2]) %calculate the T/P-Jason alias period of S2 by giving the S2 frequency (unit:1/day)
%T/P-Jason alias period: 58.74
%%
%demonstration of s_minimumLOR.m
%this example was first added in S_TIDE v1.22 update3  2021/09/03
[mL]=s_minimumLOR(1/24,{'M2';'S2';'K1';'O1';'N2';'K2';'P1';'Q1'})  
%At least 182.62 days are needed to resolve eight main constituents if sampling interval is 1/24 day (1 hour)
[mL]=s_minimumLOR(1/24,{'M2';'S2';'K1';'O1'})
%At least 14.765 days are needed to resolve four main constituents if sampling interval is 1/24 day (1 hour)
[mL]=s_minimumLOR(9.9156,{'M2';'S2';'K1';'O1';'N2';'K2';'P1';'Q1'})
%At least 9.18 years are needed to resolve eight main constituents if sampling interval is 9.9156 days (T/P-Jason satellite)

%%
%demonstration of s_quasi_HA.m
% this example was first added in S_TIDE v1.23  2021/12/13
load Qingdao.mat
 [S,H,G,coef,xout,ju,coefint,Ra,Rb]=s_quasi_HA(wl2(1:25),'1997/08/01 00:00:00',25,1,35.23,119.33);
 plot(wl2(1:25),'r.');hold on;plot(xout,'k')
%%
%demonstration of s_quasi_HA2.m
% this example was first added in S_TIDE v1.23 update  2021/12/15
 load Qingdao.mat
 for i=1:25
    time(i,:)=datevec(datenum([1997,08,01,00,00,00])+(i-1)/24);
end

[S,H,G,coef,xout,ju,Htint,Gtint,Ra,Rb]=s_quasi_HA2(wl2(1:25),time,35.23,119.33);
 %%
 %demonstration of s_equilibrium_tide.m
 % this example was first added in S_TIDE v1.23  2021/12/14
 [sumconsti,consti,amp,pha]=s_equilibrium_tide(19.5777,-156.5393,1,8760,[2005,06,28,00,00,00]);
 % this location is NDBC 51407 - Hawaii - 34 NM West of Kailua-Kona (https://www.psmsl.org/data/bottom_pressure/locations/1263.php)
 % equilibrium M2 tide  14.9cm, 46.9¡ã
 % observed    M2 tide  21.0cm, 52.8¡ã£¨https://www.psmsl.org/data/bottom_pressure/locations/data/51407_0507.tidecon£©
 
 %%
 %demonstration of s_estimate_max_tidalcurrent.m and s_tdd.m
 % this example was first added in S_TIDE v1.23 update2    2022/01/02
load tidalcurrents.mat

[Stv,Htv,Gtv,coef,xoutv,ju,Stint,Htint,Gtint,aa,bb,namev,fv,constiv]=s_tide_m8(v,1,1,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'},8,1,'spline','robustfit',20.5,[2010,08,23,01,00,00],'corrected');
[Stu,Htu,Gtu,coef,xoutu,ju,Stint,Htint,Gtint,aa,bb,nameu,fu,constiu]=s_tide_m8(u,1,1,{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'},8,1,'spline','robustfit',20.5,[2010,08,23,01,00,00],'corrected');
% 20.5 is the latitude of the observation,[2010,08,23,01,00,00] is the start time of observation, users only need to change these two inputs when perform their own harmonic analysis.
% v and u are observed tidal currents at y and x axis
[mct,angle]=s_estimate_max_tidalcurrent(Htu(:,1),Gtu(:,1),Htv(:,1),Gtv(:,1),{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'});
%mct: maximum possible tidal current speed
%angle: the angle of maximum tidal current
[tdd]=s_tdd(Htu(:,1),Gtu(:,1),{'M2';'S2';'K1';'O1';'N2';'Q1';'P1';'K2'}) %calculate the theoretical depth datum(tdd)

 %demonstration of s_rotation_spectra.m
 % this example was first added in S_TIDE v1.23 update2    2022/01/05
[px,py,pcw,pccw,freq,rdn] = s_rotation_spectra(xoutu,xoutv,1/3600);
plot(24*3600*freq,log10(pcw),'r')
hold on
plot(24*3600*freq,log10(pccw),'k')
xlim([1.8 2.2]);xlabel('cycles per day(cpd)')
ylabel('PSD(m^{2}/s^{2}/Hz)')
line([24*0.080511 24*0.080511],[-5 5],'Linestyle','--');
text(1.94,-4,'M2');legend('clockwise','counterclockwise')

%%
%demonstration of s_rtl.m
%this example was first added in S_TIDE v1.23 update3 2022/01/13
load kushiro.mat
 [rHt,rLt]=s_rtl(kushiro(1:8767)/10,[1993,01,01,00,00,00],1,42.5,6,0.1,'asymmetrical')
 %³Ë¸ß³±Ë®Î»£¨6Ð¡Ê±,ÀÛ»ýÆµÂÊ10%£©201.51cm
 %³ËµÍ³±Ë®Î»£¨6Ð¡Ê±,ÀÛ»ýÆµÂÊ10%£©148.86cm
 
 %%
 %demonstration of s_modaldecomposition
 % this example was first added in S_TIDE v1.23 update4 2022/03/29
 load data.mat; 
 % u_k1 means u speed for K1 tide
[N,Fi,dFi,mc,um]=s_modaldecomposition(depth,density,u_k1,4); % 4 baroclinic modes
%um(5,:,:) is the barotropic mode
i=1900; plot(squeeze(um(4,:,i)+um(5,:,i)+um(3,:,i)+um(2,:,i)+um(1,:,i)),depth,'k')
 set(gca,'YDir','reverse');hold on; plot(u_k1(:,i),depth,'r*')
 ylabel('Depth(m)');xlabel('u(m/s)');
 legend('reconstructed','observed')
 
 d=18;
 plot(squeeze(um(4,d,:)+um(5,d,:)+um(3,d,:)+um(2,d,:)+um(1,d,:)),'k*')
 hold on;plot(u_k1(d,:),'r');xlim([1 100])
 ylabel('u(m/s)');xlabel('Time(hour)');
 title('Depth: 85m');legend('reconstructed','observed')