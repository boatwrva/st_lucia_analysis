function [ua] = smooth_gap1( u,wt,ngap )

%revised August 23, 2023, Mark Merrifield

% Applies a low-pass filter by convolving u (input time series) with wt (filter weights).
% Gaps (NaNs) in u are filled if sum(wt(~isnan(u)))/sum(wt) > ngap.  
% The filter weights are renormalized to take into account gaps.
%

if  nargin < 2
    error(message('TooFewInputs'));
elseif nargin == 2
    ngap = 0.8;
end

u = u(:);

%set a gap flag and zero out gaps in input time series;
uflag = ones(size(u));
uflag(isnan(u)) = 0;
u_zero = u;
u_zero(isnan(u)) = 0;

%pad start and end of time series and flag series with zeros
n = length(u_zero);
nwt = length(wt);
nwt2 = floor(nwt/2);
uu = [zeros(nwt2,1);u_zero;zeros(nwt2,1)];
uuflag = [zeros(nwt2,1);uflag;zeros(nwt2,1)];

%apply convolution
ua = conv(uu,wt,'same');

%convolve flag vector with wt, and normalize output ua; 
nu = conv(uuflag,wt,'same');
ua = ua./nu;

%NaN out output if nu/sum(wt) < ngap
kn = find(nu/sum(wt) < ngap);
ua(kn) = NaN;

%trim beginning and end of output time series
ua = ua(nwt2+1:nwt2+n);

