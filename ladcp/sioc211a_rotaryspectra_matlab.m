%% rotary spectra 

%Let's load in some real data from Alford et al (2012) and have a look.
clear all

data_dir = '/home/vboatwright/OneDrive/Documents/SIO/fall23/SIO221/sioc221a/week7/'; 
fn = 'Vel_2008-2009.mat'; 
file = [data_dir fn]; 
load(file) 
%Now we will make a complex vector.
it=1:16000;
iz=min(find(Vel.z > 300));
data=Vel.u(iz,it)+i*Vel.v(iz,it);
time=Vel.dtnum(it);
dt=nanmean(diff(time));
disp(['Sample interval is ' num2str(dt*24*60) ' minutes.'])
disp(['record length is ' num2str(nanmax(time) - nanmin(time)) ' days long.'])

%% Check it for nans

disp(['The record has ' num2str(length(find(isnan(data)))) ' NaNs.'])
figure(1)
plot(time,real(data),time,imag(data))
ylabel('E/W velocity - m/s')
xlabel('Time')
legend('u','v')
datetick
%And let's zoom in:
figure(2)
plot(time,real(data),time,imag(data))
ylabel('E/W velocity - m/s')
xlabel('Time')
xlim(datenum(2009,2,1,0,0,0)+[0 5])
datetick('x','keeplimits')
legend('u','v')

%% Now we compute the rotary spectrum.
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
% amp = amp .* 2; %Now we keep the whole spectrum so we do not multiply by 2.
amp = amp / df; % this is then the definition of the spectrum
variance=nanmean(abs(data).^2)
sum_spec=sum(amp)*df
sum_spec / variance
%Check! It gives the variance.
figure(2)
loglog(f,amp,-f,amp)
legend('positive','negative')
ylabel('(m/s)^2/cpm')
xlabel('cpm')

%% windowing 

M=8; p=N/M;
datab=reshape(data,N/M,M); % this gives us an array with N/M points
% per column and M columns
b=fft(datab); % this computes the fft for each column
b=fftshift(b);
amp_b=abs(b).^2;
%% Compute frequencies - note differences!
Tb=dt*p;
dfb=1/Tb;
fn=1/2/dt;
fb=0:dfb:fn; %frequency vector, cpd, goes from 0 to Nyquist.
fb=-p/2*dfb:dfb:p/2*dfb-dfb;
%% Normalize as above
amp_b = amp_b / p.^2; % first correct for the MATLAB normalization

amp_b = amp_b / dfb; % this is then the definition of the spectrum
loglog(fb,nanmean(amp_b,2),-fb,nanmean(amp_b,2))
freqline(2*sind(50));
freqline(2);
freqline(1);
grid
legend('positive','negative')
ylabel('(m/s)^2/cpd')
xlabel('cpd')
