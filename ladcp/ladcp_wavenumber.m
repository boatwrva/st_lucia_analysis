%% SLE Wave number analysis 

data_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/processed_ship/'; 
fn = 'LADCP_station1.mat'; 

file = [data_dir fn]; 
load(file) 

%%
analysis_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/st_lucia_analysis';
addpath analysis_dir 

[nt,nz] = size(LADCP1.u);
u = LADCP1.u; v = LADCP1.v; p = LADCP1.p; 
time_dim = repmat(LADCP1.time,1,nz);

figure(1); clf; 
pcolor(time_dim,p,u); shading flat
caxis([-.25 .25]); cmocean('balance'); colorbar
shg
title('E/W velocity [m/s]')
datetick; set(gca, 'YDir', 'reverse');
xlabel('Time')
ylabel('depth [m]')

%%
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

