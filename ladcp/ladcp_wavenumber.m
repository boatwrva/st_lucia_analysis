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
clim([-.25 .25]); cmocean('balance'); colorbar
shg
title('E/W velocity [m/s]')
datetick; set(gca, 'YDir', 'reverse');
xlabel('Time')
ylabel('depth [m]')




%% make rotary spectra 

rotary = nan(nt,nz); 
column_length = [];
% make complex vector for each segment - ignore missing vals at bottom
for it = 1:nt
    iz = find(~isnan(u(it,:)));
    complex = u(it,iz)+i*v(it,iz); 
    rotary(it,iz) = complex; 
    column_length = cat(1,column_length,length(iz));
    disp(column_length)
end

% compute rotary spectra for each cast: U(t) = u(t) + i*v(t) 
% only use minimum data to preserve spectra across casts 

min_depth = min(column_length);
data = rotary(1:nt,1:min_depth);

dates = LADCP1.time; 
sampling_period = nanmean(diff(dates)); %dt between casts 
dt = minutes(sampling_period); 
N = nt; 
T=dt*N;
fundamentalf=1/T;
nyf=1/2/dt;

disp(sprintf('Sample interval is %s minutes',dt))
disp(sprintf('record length is %s long.',max(dates)-min(dates)))


dz = nanmean(diff(p,1,2),'all'); 
nz = min_depth;
depths = p(1:min_depth);
M = dz*min_depth; 
fundamentalm = 1/M; 
nym = 1/2/dz; 

disp(sprintf('Sampling interval is %d meters',dz))
disp(sprintf('record length is %d deep/long.',M))

% Check it for nans!
disp(['The record has ' num2str(length(find(isnan(data)))) ' NaNs.'])

pidx = 247; % index for depth of 997m
figure(1)
title('plotting at 1 depth: ~1000m')
plot(dates,real(data(:,pidx))); hold on
plot(dates,imag(data(:,pidx)))
ylabel('velocity [m/s]')
xlabel('Time')
legend('u','v')
datetick

castnum = 1;
figure(2)
title('plotting first cast')
subplot(1,2,1); hold on 
plot(real(data(castnum,:)),depths(castnum,:)); hold on
plot(imag(data(castnum,:)),depths(castnum,:))
set(gca,'YDir','reverse')
ylabel('depth [m]')
xlabel('velocity [m/s]')
legend('u','v')
xlim([-0.25,0.25])
subplot(1,2,2); hold on 
scatter(u(castnum,1:min_depth),depths(castnum,:)); hold on
scatter(v(castnum,1:min_depth ),depths(castnum,:))
set(gca,'YDir','reverse')
ylabel('depth [m]')
xlabel('velocity [m/s]')
legend('u','v')
xlim([-0.25,0.25])



%% Now we compute the rotary spectrum per cast .

windowing = 1; 
windows = nt; datapoints = nz;

if windowing==1 %detrend and demean
    datab = detrend(data);
end

if windowing==2 %detrend and apply a Hanning window
    w = hann(datapoints);
    w = w/sqrt(mean(w.^2));
    datab = detrend(data);
    datab = datab.*w;
end

b=fft(datab); % computes the fft for each cast
b=fftshift(b); % use fftshift to rearrange to have 0 wavenumber in the center

% you should have to make some changes for odd N 
%if rem(datapoints,2)==0
%amp_b=abs(b(:,1:datapoints/2+1)).^2; % even N
%else 
%amp_b=abs(b(:,1:(datapoints+1)/2).^2); %  odd N
%end

% square for spectrum
amp_b = abs(b(:,1:datapoints).^2); 
% normalize amp
amp_b = amp_b / datapoints.^2;  % correct for the MATLAB normalization
% amp_b = amp_b .* 2;           % keeping the whole spectrum so do not want to multiply
amp_b = amp_b / fundamentalm;   % definition of the spectrum
amp_all = mean(amp_b,1);  % average for all segemnts 

%error bars
nu=2*windows;
amp_err.upper = amp_b*nu/chi2inv(0.05/2,nu);     %upper
amp_err.lower = amp_b*nu/chi2inv(1-0.05/2,nu);   %lower 

%Check whether it gives the variance to satisfy parsevals
variance = nanmean(abs(datab).^2,'all')
sum_spec = sum(amp_all)*fundamentalm
parsevals = sum_spec/variance


% rotary wavenumber goes from negative fm to positive fm:
m=-datapoints/2*fundamentalm:fundamentalm:datapoints/2*fundamentalm-fundamentalm;

figure(3)
%loglog(m,amp_b(:,:),-m,amp_b(:,:),'LineWidth',1); hold on 
loglog(m,amp_all,-m,amp_all,'LineWidth',3)
xline(1/100,'Color','k','DisplayName','m=1/100m')
xline(1/300,'Color','k','DisplayName','m=1/300m')

legend('positive','negative')
ylabel('(m/s)^2/cpm')
xlabel('cpm')

