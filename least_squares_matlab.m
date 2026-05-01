%% least squares notes for matlab
% S Gille's notes from HW5 data2 2024 

% might need to convert the dates
date=datenum(M.data(:,1),M.data(:,2),M.data(:,3),M.data(:,4),M.data(:,5),M.data(:,6));
% store the data and convert missing values to NaN
xx=find(M.data(:,8)>9000);
data=M.data(:,8);
data(xx)=NaN;
plot(date,data,'LineWidth',2)
datetick
xlabel('date','FontSize',14)
ylabel('pressure (dbars)','FontSize',14)


% this here is how to do least squares fit in matlab 

% but we're not missing cast in station 1 so this is chill 
%tidal coefficient period (hrs)
%M2 12.42
%S2 12.00
%N2 12.66
%K1 23.93
%O1 26.87

M2=24/12.42;
S2=24/12;
N2=24/12.66;
K1=24/23.93;
O1=24/25.82;


N=length(date);
G=[ones(N,1) date-mean(date) cos(date*2*pi*M2) sin(date*2*pi*M2) ...
cos(date*2*pi*S2) sin(date*2*pi*S2) cos(date*2*pi*N2) sin(date*2*pi*N2) ...
cos(date*2*pi*K1) sin(date*2*pi*K1) cos(date*2*pi*O1) sin(date*2*pi*O1)];
xx=find(~isnan(data));
params=inv(G(xx,:)'*G(xx,:))*G(xx,:)'*data(xx)
hold on
plot(date,G*params,'LineWidth',2)
legend('data','fit','FontSize',14)

% weighted 

wt=0.2;
G2=G/wt

data2=data/wt;
uncertainties=sqrt(diag(inv(G2(xx,:)'*G2(xx,:))));
% uncertainty for mean is first value:
uncertainties(1)

standard_error_of_mean=std(data(xx))/sqrt(length(xx));
[uncertainties(1) standard_error_of_mean]

misfit = sqrt(mean((G2(xx,:)*params2-data2(xx)).^2))

