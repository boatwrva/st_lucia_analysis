%% wavenumber spectra 

% background in matlab: variance preserving spectra 

fsig=0.03; %Signal at about 30 seconds
n2s=3; % noise-to-signal ratio
dt=1; %one second
time=(1:10000)'.*dt;

% Red spectrum
b=n2s*cumsum(randn(10000,1))+cos(2*pi*fsig*time);

% White spectrum
%b=n2s*randn(10000,1)+cos(2*pi*fsig*time);
T=500;
fn=1/2/dt;
f=(0:250)./250.*fn;
bb=reshape(b,500,10000/500);
fbb=fft(detrend(bb));
amp=abs(fbb(1:251,:)).^2 / 500;
sbb=mean(amp,2);
sbb(2:250)=sbb(2:250)*2;
loglog(f,sbb,'LineWidth',3)
xlabel('Frequency (Hz)','FontSize',14)
ylabel('(units^2)/Hz','FontSize',14)

%% wavenumber spectra 

msig=0.05; %Signal at about 20 m vertical wavelength
noise_rms=0.02; % rms noise
signal = 0.05; %signal amplitude in m/s
dz=1/2; %1/2 m spacing - 5000 m profile
z=(1:10000)'.*dz;

% White spectrum
b=noise_rms*randn(10000,1)+signal.*cos(2*pi*msig*z);
L=250;
mn=1/2/dz; %Nyquist wavenumber
m=(0:250)./250.*mn; %wavenumber vector
bb=reshape(b,500,10000/500);
fbb=fft(detrend(bb));
amp=abs(fbb(1:251,:)).^2 / 500;
sbb=mean(amp,2);
sbb(2:250)=sbb(2:250)*2;
subplot(121)
plot(bb(:,1),z(1:500));
axis ij
xlabel('m/s')
ylabel('depth / m')
subplot(122)
loglog(m,sbb,'LineWidth',3)
xlabel('Wavenumber (cycles per meter)','FontSize',14)
ylabel('(m/s)^2/cpm','FontSize',14)



%% wavenumber-frequency spectra 

H=1000;
dt=0.05; %time interval about 20 times a day
dz=10; %10-m
t=dt:dt:10; %time in days
z=(dz:dz:H)'; %depth in m
%indices for the calculation
iz=1:length(z);
it=1:length(t);
%Generate a propagating signal with no noise.
[tt,zz]=meshgrid(t,z);
sig1=sin(3*pi*zz/H + 2*pi*1/3*tt);
%Change the sign to make propagation go the other way.
data=sig1;
% plot
figure(2)
imagesc(t,z,data);colorbar; axis ij
xlabel('time (days)')
ylabel('depth (m)')

% compute frequency-wavenumber spectra 
% 1) fourier transform in 2D 

[m,n]=size(data);
fn=1/2/dt;
kn=1/2/dz;

%fundamental frequency and wavenumber
df=1./n./dt;
dk=1./m./dz;
% make frequencies and wavenumbers that run from -Nyquist to + Nyquist
f=[-fliplr(1:(n/2)) 0 (1:(n/2-1))].*df;
k=[-fliplr(1:(m/2)) 0 (1:(m/2-1))].'.*dk;

% Fourier transform in two dimensions
% here we use fft2 for the 2-d Fourier transform
% and fftshift to reorder the Fourier transform
st=fftshift(fft2(data))/m/n;

% turn this into a spectrum
spec=st.*conj(st)./df./dk; %UNITS: (m/s)^2/cpd/cpm

% and plot
figure(3)
imagesc(f,k,log10(spec)); axis xy
colormap(jet)
shg
xlabel('\omega / cpd')
ylabel('m / cpm')

% this plots the full 2-d plane
% But you might want a half plane. In that case, you should scale by a factor of 2:
figure(4)
spec=2*spec(:,101:200);
imagesc(f(101:200),k,log10(spec)); axis xy
colormap(jet)
shg
xlabel('\omega / cpd')
ylabel('m / cpm')


%% segmenting in space - time 

t=dt:dt:30; %time in days
z=(dz:dz:H)'; %depth in m
%indices for the calculation
iz=1:length(z);
it=1:length(t);
%Generate a propagating signal with some noise.
[tt,zz]=meshgrid(t,z);
sig1=sin(3*pi*zz/H - 2*pi*1/3*tt) + .5*randn(length(z),length(t));
data=sig1;
% plot
figure(2)
imagesc(t,z,data);colorbar; axis ij
xlabel('time (days)')
ylabel('depth (m)')

% segmenting 

icount=0; clear st;
for i=1:50:501
icount=icount+1;
data_use=data(:,i:i+99);
n_use=size(data_use,2);
st(:,:,icount)=fftshift(fft2(data_use))/m/n_use;
end
spec=sum(abs(st).^2/df/dk,3);
f=[-fliplr(1:(n_use/2)) 0 (1:(n_use/2-1))].*df;
figure(3)
imagesc(f,k,log10(spec)); axis xy
colormap(jet)
shg
xlabel('\omega / cpd')
ylabel('m / cpm')



%% using real data 

data_dir = '/home/vboatwright/OneDrive/Documents/SIO/fall23/SIO221/sioc221a/week7/'; 
fn = 'Vel_2008-2009.mat'; 

file = [data_dir fn]; 
load(file) 


figure(1); clf; 
pcolor(Vel.dtnum,Vel.z,Vel.u); shading flat
caxis([-.25 .25]);colormap(turbo);colorbar
shg
title('E/W velocity - m/s')
datetick
xlabel('Time')
ylabel('depth / m')

figure(2)
iz=min(find(Vel.z > 300));
it=1:length(Vel.dtnum);
data=Vel.u(iz,it);
time=Vel.dtnum(it);
plot(time,data)
ylabel('E/W velocity - m/s')
xlabel('Time')
datetick


%Now we will make a complex vector.
it=1:16000;
data=Vel.u(iz,it)+i*Vel.v(iz,it);
time=Vel.dtnum(it);
dt=nanmean(diff(time));
disp(['Sample interval is ' num2str(dt*24*60) ' minutes.'])
disp(['record length is ' num2str(nanmax(time) - nanmin(time)) ' days long.'])


% Check it for nans!
disp(['The record has ' num2str(length(find(isnan(data)))) ' NaNs.'])

%Now we compute the rotary spectrum.
data=detrend(data);
a=fft(data); %Note we keep all the coefficients instead of throwing out half.
a=fftshift(a); %And now we use fftshift to rearrange them to have 0 frequency in the middle.
N=length(data);
amp=abs(a).^2; % for even N

%% which freqs do these correspond to?
T=dt*N;
df=1/T;
fn=1/2/dt;
%Now frequency goes from negative fn to positive fn:
f=-N/2*df:df:N/2*df-df;

amp = amp / N.^2; % first correct for the MATLAB normalization
%amp = amp .* 2; %Now we keep the whole spectrum so we do not multiply by 2.
amp = amp / df; % this is then the definition of the spectrum
variance=nanmean(abs(data).^2)
sum_spec=sum(amp)*df
sum_spec / variance
%Check! It gives the variance.
loglog(f,amp,-f,amp)
legend('positive','negative')
ylabel('(m/s)^2/cpd')
xlabel('cpd')

%% window for precision 

M=8; p=N/M;
datab=reshape(data,N/M,M); % this gives us an array with N/M points
% per column and M columns
b=fft(datab); % this computes the fft for each column
b=fftshift(b);
amp_b=abs(b).^2;

% Compute frequencies - note differences!
Tb=dt*p;
dfb=1/Tb;
fn=1/2/dt;
fb=0:dfb:fn; %frequency vector, cpd, goes from 0 to Nyquist.
fb=-p/2*dfb:dfb:p/2*dfb-dfb;

% Normalize as above
amp_b = amp_b / p.^2; % first correct for the MATLAB normalization
amp_b = amp_b / dfb; % this is then the definition of the spectrum
loglog(fb,nanmean(amp_b,2),-fb,nanmean(amp_b,2))
freqline(2*sind(50));
freqline(2);
freqline(1);
grid()
legend('positive','negative')
ylabel('(m/s)^2/cpd')
xlabel('cpd')
