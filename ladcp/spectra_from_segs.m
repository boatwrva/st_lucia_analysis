function [fb,amp_b,amp_err] = spectra_from_segs(datab,dt,win)
%inputs: 
%datab: (segmented data)
%dt: sample period
%win: hanning or detrend

%outputs: 
%fb: freqeuncy vector 
%amp_b: amplitudes
%amo_err: upper and lower bounds on amp error


[p,M]=size(datab);

if win==1 %detrend and demean
    datab = detrend(datab);
end

if win==2 %detrend and apply a Hanning window
    w = hann(p);
    w = w/sqrt(mean(w.^2));
    datab = detrend(datab);
    datab = datab.*w;
end

b=fft(datab); % computes the fft for each column
if rem(p,2)==0
amp_b=abs(b(1:p/2+1,:)).^2; % even N
else 
amp_b=abs(b(1:(p+1)/2,:).^2); %  odd N
end

% frequency
Tb=dt*p;     % Length
dfb=1/Tb;    % Fundamental
fn=1/2/dt;   % Nyquist 
fb=0:dfb:fn; % frequency vector

% normalize amp
amp_b = amp_b / p.^2;   % correct for the MATLAB normalization
amp_b = amp_b .* 2;     % correct for the lost variance.
amp_b = amp_b / dfb;    % definition of the spectrum
amp_b = mean(amp_b,2);  % average for all segemnts 

%error bars
nu=2*M;
amp_err.upper = amp_b*nu/chi2inv(0.05/2,nu);     %upper
amp_err.lower = amp_b*nu/chi2inv(1-0.05/2,nu);   %lower 


% parsevals 
% variance=mean(data.^2);
% sum_spec=mean(sum(amp_b)*dfb);
% parsevals_b = sum_spec / variance;
end