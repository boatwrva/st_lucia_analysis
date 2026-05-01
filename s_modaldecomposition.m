function [N,Fi,dFi,mc,um]=s_modaldecomposition(depth,density,u,m)
%written by Haidong Pan,First Institute of Oceanography
%2022/03/29,panhaidong_phd@qq.com
%Dynamical Modal Decomposition for internal tides
%input: depth:  depth of current records (unit:m) 
%      density: density of ocean water at given depth (kg/m^3)
%      u: current records at given depth of one constituent (missing values labeled as NaN)
%      m: number of baroclinic modes
%output: N: Buoyancy frequency (/s)
%       Fi:eigenfunction
%       dFi: Derivative of eigenfunction
%       mc: modal coefficients
%       um: tidal currents of each mode (m+1 modes, m baroclinic modes, while the last mode is barotropic mode)
%Citing:
% Pan et al.(2018)  Exploration of tidal-fluvial interaction in the Columbia River estuary using S_TIDE
% Cao et al.(2015) Extraction of Internal Tidal Currents and Reconstruction of Full-Depth Tidal Currents from Mooring Observations
% see https://www.bilibili.com/video/BV15L411K7iG?spm_id_from=333.880 for more details
% see s_demo.m for examples
disp(['If you have any ideas or questions about S_TIDE,please send emails to me, I am happy to discuss with you'])
disp(['********** Dynamical Modal Decomposition****************'])
disp(['************* written by Haidong Pan(FIO) ********************'])
mh=max(depth);di=density;
N=sqrt(diff(di)./diff(depth)*9.81./di(1:end-1));%Calculate Buoyancy frequency¼ÆËã¸¡ÐÔÆµÂÊ
N(end+1)=N(end);
dh=diff(depth);dh(end+1)=dh(end);
%calculate eigenfunctions (WKB approximation used) ¼ÆËã´¹Ïò±¾Õ÷º¯Êý
for j=1:m
for i=1:length(depth)
    Nz= sum(N(1:i).*dh(1:i));%dh*sum(N(1:i));
    Fi(i,j)=(-1)^j/j*sqrt(nanmean(N)/N(i))*sin(pi*j/nanmean(N)/mh*Nz);
    dFi(i,j)=(-1)^j/mh*pi*sqrt(N(i)/nanmean(N))*cos(pi*j/nanmean(N)/mh*Nz);
end
end
dFi(:,m+1)=1/mh;%ÕýÑ¹Ä£Ì¬  barotropic mode
for i=1:length(u(1,:))
  Y=u(:,i);gd=find(isfinite(Y));
  K=dFi(gd,:);
  mc(:,i)=inv(K'*K)*K'*Y(gd);%¼ÆËãÄ£Ì¬ÏµÊýcalculate modal coefficients 
end
 for i=1:m+1
     for j=1:length(depth)
         for k=1:length(u(1,:))
            um(i,j,k)=mc(i,k)*dFi(j,i); %¹¹½¨¸÷Ä£Ì¬Á÷ËÙ consturct tidal currents of each mode
         end
     end
 end
    