% file for processing initial LADCP from Santa Lucia 

dir_in = '/media/vboatwright/KBZ/SR2503_cruise/';
dir_ladcp = '/media/vboatwright/KBZ/SR2503_scienceparty_share/LADCP/processed/Survey_Stations/';
stn = 11; 
file = sprintf('%03d.mat',stn);
fn = [dir_ladcp file]; 
S = load(fn); 

%% check TS 

temp = S.dr.ctd_t; 
sal = S.dr.ctd_s; 

figure()
plot(temp)
figure()
plot(sal)

%% 

casts = 1:12; 
numcasts = length(casts);
stations = [1,5,9,11,7,2,6,10,12,8,4,3]; 
deepest = 1200; 
u_velocity = nan(12,deepest); 
v_velocity = nan(12,deepest); 
z_depth = nan(12,deepest); 

figure(1); hold on 
ii = 0;
for icast=2:13
    ii = ii + 1; 
    file = sprintf('%03d.mat',icast);
    fn = [dir_ladcp file]; 
    S = load(fn);
    data = S.dr; 
    zlevs = length(data.z);
    u_velocity(ii,1:zlevs) = data.u; 
    v_velocity(ii,1:zlevs) = data.v; 
    z_depth(ii,1:zlevs) = data.z; 

    subplot(3,4,ii); hold on 
    x = zeros(zlevs,1);
    y = zeros(zlevs,1);
    w = zeros(zlevs,1); 

    title(sprintf('Station %d',stations(icast-1)))
    quiver3(x,y,data.z,data.u,data.v,w,0)
    zlabel('Depth [m]'); set(gca, 'ZDir', 'reverse');
    xlabel('Velocity [m/s]')

end


%% save all casts 

ladcp_casts = struct; 
ladcp_casts.u_vel = u_velocity; 
ladcp_casts.v_vel = v_velocity; 
ladcp_casts.z = z_depth; 
ladcp_casts.stations = stations;

%% organize by location in grid 

castnumbers = [1:13]; 
stationnumbers = [0,1,5,9,11,7,2,6,10,12,8,4,3]; 
% going from 1 at northwest ---> 4 at northeast; 5 at middle west ---> ; 12 at southeast 
cardinallocs = [0,9,5,1,3,7,10,6,2,4,8,12,11]; 
stationlocs = cardinallocs(2:end);

% plot with u,v both showing (not 3D) 
figure(); hold on; set(gcf,'Color','White') 
ii = 0;
for icast=1:12
    subplot(3,4,stationlocs(icast)); hold on 
    title(sprintf('Station %d',stations(icast)))
    
    pu = plot(u_velocity(icast,:),z_depth(icast,:),'Color','#A2142F','DisplayName','u'); hold on 
    pv = plot(v_velocity(icast,:),z_depth(icast,:),'Color','#0072BD','DisplayName','v'); hold on 
    xline(0,'Color','Red','LineStyle','--','DisplayName','Zero')
    % z in depth (0 is surface) 
    xlabel('Velocity [m/s]')
    ylabel('Depth [m]'); set(gca,'YDir','reverse')
    legend()

end


   
%% again in 3d 

figure(1); hold on 
ii = 0;
zlevs = 1200;
for icast=1:12
    subplot(3,4,stationlocs(icast)); hold on 
    title(sprintf('Station %d',stations(icast)))
    
    x = zeros(1,zlevs);
    y = zeros(1,zlevs);
    w = zeros(1,zlevs); 

    quiver3(x,y,z_depth(icast,:),u_velocity(icast,:),v_velocity(icast,:),w,0)
    zlabel('Depth [m]'); set(gca, 'ZDir', 'reverse');
    xlabel('U Velocity [m/s]')
    ylabel('V Velocity [m/s]')

end



%% quiver plot for all stations 


figure(1); hold on 
for pp=1:6
    subplot(3,4,pp); hold on 
    x = zeros(length(z_uv),1);
    y = zeros(length(z_uv),1);
    w = zeros(length(z_uv),1); 

    u = u_velocity(pp,:); 
    v = v_velocity(pp,:); 
    z = z_depth(pp,:); 
    title(sprintf('Station %d',stations(pp)))
    quiver3(x,y,z,u,v,w,0)
end

%%
data = S.dr; 

% current velocities of up/down profiles

z_uv = data.z; 

u_down = data.u_do; 
v_down = data.v_do; 

u_up = data.u_up; 
v_up = data.v_up; 

u = data.u; 
v = data.v; 

u_bot = data.ubot; 
v_bot = data.vbot; 
z_bot = data.zbot; 

adcp_time = data.tim; 
adcp_hours = data.tim_hour; 

lons = data.shiplon; 
lats = data.shiplat; 

%dnum = data.datenum; 
%dates = datetime(dnum,'ConvertFrom','datenum');

%% quiver plot of one station to compare u/v, other data 

x = zeros(length(z_uv),1);
y = zeros(length(z_uv),1);
w = zeros(length(z_uv),1); 


figure(1)
subplot(3,4,1)
title(sprintf('Station %d',stn))
quiver3(x,y,z_uv,u,v,w,0)

subplot(3,4,2)
title(sprintf('Station %d',stn))
quiver3(x,y,z_uv,u,v,w,0)

subplot(3,4,3)
title(sprintf('Station %d',stn))
quiver3(x,y,z_uv,u,v,w,0)

subplot(3,4,4)
title(sprintf('Station %d',stn))
quiver3(x,y,z_uv,u,v,w,0)

%%

t1 = 200; 
subplot(1,3,3)
title(sprintf('time %s',dates(t1)))
nd = datenum(dates(t1:t1+5));
u_avg = mean(u(:,t1:t1+10),2);
v_avg = mean(v(:,t1:t1+10),2); 
subplot(1,2,1)
quiver3(x,y,z_uv,u(:,t1),v(:,t1),w,0)
subplot(1,2,2)
quiver3(x,y,z_uv,u_avg,v_avg,w,0)



%% 3D plot 


f = figure(1); clf; set(gcf,'Color','w'); hold on 
ax(1) = subplot(2,1,1); hold on 

p1 = pcolor(adcp_time,z_uv,u); %,'DisplayName','Zonal Velocity'); %'Color',orange);
xlabel('Time'); ylabel('Depth [m]');
shading flat;
set(gca, 'YDir', 'reverse');
c = colorbar;
c.Label.String = 'Velocity [m/s]';
datetick('x','mm/yy','keepticks','keeplimits')
title('Zonal Velocity')


ax(2) = subplot(2,1,2); hold on 
p2 = pcolor(dates,z_uv,v); %,'DisplayName','Meridional Velocity'); %,'Color',orange);
shading flat
set(gca, 'YDir', 'reverse');
xlabel('Time'); ylabel('Depth [m]');
c = colorbar;
c.Label.String = 'Velocity [m/s]';
datetick('x','mm/yy','keepticks','keeplimits')
title('Meridional Velocity')


%% problem 1: pick a depth and look at how the u and v components of velocity rotate in time
% is this consistent with what we discussed in class for internal waves? 
% if you like (and know how) you can also bandpass it around a frequency of your choice to see that more easily.

% depth values start at 130m depth, index 70

u500 = u(z_uv==500,:); 
v500 = v(z_uv==500,:); 

u1000 = u(z_uv==1000,:); 
v1000 = v(z_uv==1000,:); 

u3000 = u(z_uv==3000,:); 
v3000 = v(z_uv==3000,:); 



f = figure(2); clf; set(gcf,'Color','w'); hold on 
ax(1) = subplot(3,1,1); hold on 

p1 = plot(dates,u500,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates,v500,'DisplayName','Meridional Velocity'); %'Color',orange);
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
datetick('x','mm/yy','keepticks','keeplimits')
title('Velocity at 500m depth')


ax(2) = subplot(3,1,2); hold on 

p1 = plot(dates,u1000,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates,v1000,'DisplayName','Meridional Velocity'); %'Color',orange);
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
datetick('x','mm/yy','keepticks','keeplimits')
title('Velocity at 1000m depth')

ax(3) = subplot(3,1,3); hold on 

p1 = plot(dates,u3000,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates,v3000,'DisplayName','Meridional Velocity'); %'Color',orange);
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
datetick('x','mm/yy','keepticks','keeplimits')
title('Velocity at 3000m depth')


idx0 = 1; idx1 = 100;


f = figure(3); clf; set(gcf,'Color','w'); hold on 
ax(1) = subplot(3,1,1); hold on 

p1 = plot(dates(idx0:idx1),u500(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v500(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
title('Velocity at 500m depth')


ax(2) = subplot(3,1,2); hold on 

p1 = plot(dates(idx0:idx1),u1000(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v1000(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
title('Velocity at 1000m depth')

ax(3) = subplot(3,1,3); hold on 

p1 = plot(dates(idx0:idx1),u3000(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v3000(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
title('Velocity at 3000m depth')



%% problem 2: multiple times and multiple depths 


u600 = u(z_uv==600,:); 
v600 = v(z_uv==600,:); 

u800 = u(z_uv==800,:); 
v800 = v(z_uv==800,:); 

idx0 = 1; idx1 = 100;

f = figure(1); clf; set(gcf,'Color','w'); hold on 
ax(1) = subplot(4,1,1); hold on 

p1 = plot(dates(idx0:idx1),u500(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v500(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
%xtickformat('x','mm/yy HH:MM','keepticks','keeplimits')
title('Velocity at 500m depth')


ax(2) = subplot(4,1,2); hold on 

p1 = plot(dates(idx0:idx1),u600(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v600(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
title('Velocity at 600m depth')

ax(3) = subplot(4,1,3); hold on 

p1 = plot(dates(idx0:idx1),u800(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v800(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
title('Velocity at 800m depth')

ax(4) = subplot(4,1,4); hold on 

p1 = plot(dates(idx0:idx1),u1000(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v1000(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
title('Velocity at 1000m depth')




idx0 = 100; idx1 = 200;

f = figure(5); clf; set(gcf,'Color','w'); hold on 
ax(1) = subplot(4,1,1); hold on 

p1 = plot(dates(idx0:idx1),u500(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v500(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
title('Velocity at 500m depth')


ax(2) = subplot(4,1,2); hold on 

p1 = plot(dates(idx0:idx1),u600(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v600(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
title('Velocity at 600m depth')

ax(3) = subplot(4,1,3); hold on 

p1 = plot(dates(idx0:idx1),u800(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v800(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
title('Velocity at 800m depth')

ax(4) = subplot(4,1,4); hold on 

p1 = plot(dates(idx0:idx1),u1000(idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(dates(idx0:idx1),v1000(idx0:idx1),'DisplayName','Meridional Velocity'); %'Color',orange);
yline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()
title('Velocity at 1000m depth')




%% 3d plot to see depth - but zoom in time 

idx0 = 1; idx1 = 100; 

f = figure(1); clf; set(gcf,'Color','w'); hold on 
ax(1) = subplot(2,1,1); hold on 

p1 = pcolor(dates(idx0:idx1),z_uv,u(:,idx0:idx1)); %,'DisplayName','Zonal Velocity'); %'Color',orange);
xlabel('Time'); ylabel('Depth [m]');
shading flat; grid on 
set(gca, 'GridLineStyle', '--','GridAlpha', 0.7,'Layer','top') 
set(gca, 'YDir', 'reverse');
c = colorbar;
c.Label.String = 'Velocity [m/s]';
datetick('x','mm/yy','keepticks','keeplimits')
title('Zonal Velocity')


ax(2) = subplot(2,1,2); hold on 
p2 = pcolor(dates(idx0:idx1),z_uv,v(:,idx0:idx1)); %,'DisplayName','Meridional Velocity'); %,'Color',orange);
shading flat; grid on 
set(gca, 'GridLineStyle', '--','GridAlpha', 0.7,'Layer','top'); % Slight transparency

set(gca, 'YDir', 'reverse');
xlabel('Time'); ylabel('Depth [m]');
c = colorbar;
c.Label.String = 'Velocity [m/s]';
datetick('x','mm/yy','keepticks','keeplimits')
title('Meridional Velocity')

%%
% lets try to put them on the same plot 
num_dates = datenum(dates(idx0:idx1));

f = figure(6); clf; set(gcf,'Color','w'); hold on 
p1 = pcolor(num_dates,z_uv,u(:,idx0:idx1),'DisplayName','Zonal Velocity'); %'Color',orange);
shading flat; set(gca, 'YDir', 'reverse'); 
c = colorbar;
c.Label.String = 'Velocity [m/s]';
hold on 
[C,h] = contour(num_dates,z_uv,v(:,idx0:idx1),15,'DisplayName','Meridional Velocity','Linecolor','k'); 
%clabel(C, h, 'FontSize', 10, 'Color', 'red'); % Customize label appearance

set(gca, 'YDir', 'reverse'); hold on 
xlabel('Time'); ylabel('Depth [m]');

datetick('x','mm/yy','keepticks','keeplimits')
title('Subsurface Velocity')
legend({'Zonal Velocity', 'Meridional Velocity'}, 'Location', 'northeast');


%% oops: plot a few timesteps 

t1 = 100; 
f = figure(6); clf; set(gcf,'Color','w'); hold on 
title('Velocities Across Depth')
subplot(1,3,1); hold on
subtitle(sprintf('Time @ %s',dates(t1)))
p1 = plot(u(:,t1), z_uv,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(v(:,t1), z_uv,'DisplayName','Meridional Velocity'); %'Color',orange);
%set(gca, 'YDir', 'reverse'); 
xline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()

t1 = 200;
subplot(1,3,2); hold on 
subtitle(sprintf('Time @ %s',dates(t1)))
p1 = plot(u(:,t1), z_uv,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(v(:,t1), z_uv,'DisplayName','Meridional Velocity'); %'Color',orange);
%set(gca, 'YDir', 'reverse'); 
xline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()

t1 = 300;
subplot(1,3,3); hold on 
subtitle(sprintf('Time @ %s',dates(t1)))
p1 = plot(u(:,t1), z_uv,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(v(:,t1), z_uv,'DisplayName','Meridional Velocity'); %'Color',orange);
%set(gca, 'YDir', 'reverse'); 
xline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()


%% quiver plot because victoria small brain 

t1 = 100; 
f = figure(6); clf; set(gcf,'Color','w'); hold on 
title('Velocities Across Depth')
subplot(1,3,1); hold on
subtitle(sprintf('Time @ %s',dates(t1)))
p1 = plot(u(:,t1), z_uv,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(v(:,t1), z_uv,'DisplayName','Meridional Velocity'); %'Color',orange);
set(gca, 'YDir', 'reverse'); 
xline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()

t1 = 200;
subplot(1,3,2); hold on 
subtitle(sprintf('Time @ %s',dates(t1)))
p1 = plot(u(:,t1), z_uv,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(v(:,t1), z_uv,'DisplayName','Meridional Velocity'); %'Color',orange);
set(gca, 'YDir', 'reverse'); 
xline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()

t1 = 300;
subplot(1,3,3); hold on 
subtitle(sprintf('Time @ %s',dates(t1)))
p1 = plot(u(:,t1), z_uv,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(v(:,t1), z_uv,'DisplayName','Meridional Velocity'); %'Color',orange);
set(gca, 'YDir', 'reverse'); 
xline(0,'Color','red','Linestyle','--','DisplayName','Zero')
xlabel('Time'); ylabel('Velocity [m/s]'); legend()

% quiver3(Z,U,V,W)
%quiver3(z_uv,u(:,t1),v(:,t1))



%% 
x = zeros(length(z_uv),1);
y = zeros(length(z_uv),1);
w = zeros(length(z_uv),1); 


figure()
subplot(1,3,1);  
t1 = 100; 
nd = datenum(dates(t1:t1+5));
%quiver(nd,z_uv,u(:,t1:t1+5),v(:,t1:t1+5),0)

quiver3(x,y,z_uv,u(:,t1),v(:,t1),w,0)
set(gca, 'ZDir', 'reverse'); 
xlabel('u [m/s]'); ylabel('v [m/s]'); zlabel('z [m]'); 
title(sprintf('t = %s',dates(t1)))



u_avg = mean(u(:,t1:t1+10),2);
v_avg = mean(v(:,t1:t1+10),2); 
%quiver3(x,y,z_uv,u_avg,v_avg,w,0)

subplot(1,3,2); 
t1 = 150; 
nd = datenum(dates(t1:t1+5));
quiver3(x,y,z_uv,u(:,t1),v(:,t1),w,0)
xlabel('u [m/s]'); ylabel('v [m/s]'); zlabel('z [m]')
set(gca, 'ZDir', 'reverse'); 
title(sprintf('t = %s',dates(t1)))

t1 = 200; 
subplot(1,3,3); 
% nd = datenum(dates(t1:t1+5));
quiver3(x,y,z_uv,u(:,t1),v(:,t1),w,0)
set(gca, 'ZDir', 'reverse'); 
xlabel('u [m/s]'); ylabel('v [m/s]'); zlabel('z [m]'); 
title(sprintf('t = %s', datestr(dates(t1)))); % Convert datetime to string


%%

figure()
t1 = 100; 
title(sprintf('time %s',dates(t1)))
nd = datenum(dates(t1:t1+5));
%quiver(nd,z_uv,u(:,t1:t1+5),v(:,t1:t1+5),0)
u_avg = mean(u(:,t1:t1+10),2);
v_avg = mean(v(:,t1:t1+10),2); 
subplot(1,2,1)
quiver3(x,y,z_uv,u(:,t1),v(:,t1),w,0)
subplot(1,2,2)
quiver3(x,y,z_uv,u_avg,v_avg,w,0)

figure()
t1 = 150; 
title(sprintf('time %s',dates(t1)))
nd = datenum(dates(t1:t1+5));
u_avg = mean(u(:,t1:t1+10),2);
v_avg = mean(v(:,t1:t1+10),2); 
subplot(1,2,1)
quiver3(x,y,z_uv,u(:,t1),v(:,t1),w,0)
subplot(1,2,2)
quiver3(x,y,z_uv,u_avg,v_avg,w,0)



t1 = 200; 
subplot(1,3,3)
title(sprintf('time %s',dates(t1)))
nd = datenum(dates(t1:t1+5));
u_avg = mean(u(:,t1:t1+10),2);
v_avg = mean(v(:,t1:t1+10),2); 
subplot(1,2,1)
quiver3(x,y,z_uv,u(:,t1),v(:,t1),w,0)
subplot(1,2,2)
quiver3(x,y,z_uv,u_avg,v_avg,w,0)


%% stretch coordinates 

% N_m z_m = integral (N(z))dz 
% find N from N2 

N = sqrt(n2);

figure(); 
plot(N,z_n2)
set(gca, 'YDir', 'reverse'); 


% discrete integral N2(z) 

intN = nan(size(N));
integralN = 0; 

for ii = 1:length(N)-1
    intN(ii) = (N(ii+1)-N(ii)) * (z_n2(ii+1) - z_n2(ii)); 
    integralN = integralN + intN(ii);
end

disp(integralN)

figure(); hold on
plot(N,z_n2,'DisplayName','N')
plot(intN,z_n2,'DisplayName','Int(N)')
xline(0,'Color','red','Linestyle','--')
set(gca, 'YDir', 'reverse'); legend()

%%


% calculate shear 
[nz,nt] = size(u);
shear_u = nan(nz,nt);
shear_v = nan(nz,nt);

for ii = 1:nz-1
    shear_u(ii,:) = (u(ii+1,:)-u(ii,:)) / (z_uv(ii+1) - z_uv(ii)); 
    shear_v(ii,:) = (v(ii+1,:)-v(ii,:)) / (z_uv(ii+1) - z_uv(ii)); 
end


%% stretching coordinates 

% reverse to go from bottom to top when cumulative summing
csumN = cumsum(N,'omitnan','reverse');

% take real only 
rsumN = real(csumN);

n_m = max(N,[],'omitnan');
z_m = rsumN/n_m; 

% z_m = intN2/ n2_mean; 


figure()
plot(z_m,z_n2)
title('Stretched Coordinates')
xlabel('Stretched z'); ylabel('Z [m]')
set(gca, 'YDir', 'reverse'); 


%% now lets plot in stretched coordinate space 

u_interp = interp1(z_uv,u,z_n2);
v_interp = interp1(z_uv,v,z_n2);


idx0 = 1; idx1 = 100; 

f = figure(); clf; set(gcf,'Color','w'); hold on 
ax(3) = subplot(2,2,3); hold on 

p1 = pcolor(dates(idx0:idx1),z_m,u_interp(:,idx0:idx1)); %,'DisplayName','Zonal Velocity'); %'Color',orange);
xlabel('Time'); ylabel('Stretched Depth [z_m = int Ndz / N_m]');
shading flat; grid on; set(gca, 'GridLineStyle', '--','GridAlpha', 0.7,'Layer','top') 
%set(gca, 'YDir', 'reverse');
c = colorbar; c.Label.String = 'Velocity [m/s]';
title('Zonal Velocity')


ax(4) = subplot(2,2,4); hold on 
p2 = pcolor(dates(idx0:idx1),z_m,v_interp(:,idx0:idx1)); %,'DisplayName','Meridional Velocity'); %,'Color',orange);
xlabel('Time'); ylabel('Stretched Depth [z_m = int Ndz / N_m]');
shading flat; grid on; set(gca, 'GridLineStyle', '--','GridAlpha', 0.7,'Layer','top'); % Slight transparency
%set(gca, 'YDir', 'reverse');
c = colorbar; c.Label.String = 'Velocity [m/s]';
title('Meridional Velocity')



% compare with real coordinates 

ax(1) = subplot(2,2,1); hold on 

p1 = pcolor(dates(idx0:idx1),z_n2,u_interp(:,idx0:idx1)); %,'DisplayName','Zonal Velocity'); %'Color',orange);
xlabel('Time'); ylabel('Depth [m]');
shading flat; grid on; set(gca, 'GridLineStyle', '--','GridAlpha', 0.7,'Layer','top') 
set(gca, 'YDir', 'reverse');
c = colorbar; c.Label.String = 'Velocity [m/s]';
title('Zonal Velocity')


ax(2) = subplot(2,2,2); hold on 
p2 = pcolor(dates(idx0:idx1),z_n2,v_interp(:,idx0:idx1)); %,'DisplayName','Meridional Velocity'); %,'Color',orange);
shading flat; grid on; set(gca, 'GridLineStyle', '--','GridAlpha', 0.7,'Layer','top'); % Slight transparency
set(gca, 'YDir', 'reverse');
xlabel('Time'); ylabel('Depth [m]');
c = colorbar; c.Label.String = 'Velocity [m/s]';
title('Meridional Velocity')


%% look at a few timesteps across stretched depth coordinate 

set(groot,'defaultAxesFontSize',18)
%fontsize(16,"default")

t1 = 100; 
f = figure(6); clf; set(gcf,'Color','w'); hold on 
title('Velocities Across Stretched Depth')
subplot(1,3,1); hold on
subtitle(sprintf('Time @ %s',dates(t1)))
p1 = plot(u_interp(:,t1), z_m,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(v_interp(:,t1), z_m,'DisplayName','Meridional Velocity'); %'Color',orange);
%set(gca, 'YDir', 'reverse'); 
xline(0,'Color','red','Linestyle','--','DisplayName','Zero')
ylabel('Stretched z_m'); xlabel('Velocity [m/s]'); legend()

t1 = 200;
subplot(1,3,2); hold on 
subtitle(sprintf('Time @ %s',dates(t1)))
p1 = plot(u_interp(:,t1), z_m,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(v_interp(:,t1), z_m,'DisplayName','Meridional Velocity'); %'Color',orange);
%set(gca, 'YDir', 'reverse'); 
xline(0,'Color','red','Linestyle','--','DisplayName','Zero')
ylabel('Stretched z_m'); xlabel('Velocity [m/s]'); legend()

t1 = 300;
subplot(1,3,3); hold on 
subtitle(sprintf('Time @ %s',dates(t1)))
p1 = plot(u_interp(:,t1), z_m,'DisplayName','Zonal Velocity'); %'Color',orange);
p2 = plot(v_interp(:,t1), z_m,'DisplayName','Meridional Velocity'); %'Color',orange);
%set(gca, 'YDir', 'reverse'); 
xline(0,'Color','red','Linestyle','--','DisplayName','Zero')
ylabel('Stretched z_m'); xlabel('Velocity [m/s]'); legend()


%% shear in stretched coordinate 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dudz_interp = interp1(z_uv,shear_u,z_n2);
dvdz_interp = interp1(z_uv,shear_v,z_n2);

%%

% calculate shear 
[nzm,nt] = size(u_interp);
dudz_interp = nan(size(u_interp));
dudz_interp = nan(size(v_interp));

for ii = 1:nzm-1
    dudz_interp(ii,:) = (u_interp(ii+1,:)-u_interp(ii,:)) / (z_m(ii+1) - z_m(ii)); 
    dvdz_interp(ii,:) = (v_interp(ii+1,:)-v_interp(ii,:)) / (z_m(ii+1) - z_m(ii)); 
end

%%

idx0 = 100; idx1 = 600;

f = figure(); clf; set(gcf,'Color','w'); hold on 
ax(3) = subplot(2,2,3); hold on 

p1 = pcolor(dates(idx0:idx1),z_m,dudz_interp(:,idx0:idx1)); %,'DisplayName','Zonal Velocity'); %'Color',orange);
xlabel('Time'); ylabel('Stretched Depth [z_m = int Ndz / N_m]');
shading flat; grid on; set(gca, 'GridLineStyle', '--','GridAlpha', 0.7,'Layer','top') 
%set(gca, 'YDir', 'reverse');
c = colorbar; c.Label.String = 'Shear [1/s]';
title('du/dz')


ax(4) = subplot(2,2,4); hold on 
p2 = pcolor(dates(idx0:idx1),z_m,dvdz_interp(:,idx0:idx1)); %,'DisplayName','Meridional Velocity'); %,'Color',orange);
shading flat; grid on; set(gca, 'GridLineStyle', '--','GridAlpha', 0.7,'Layer','top'); % Slight transparency
xlabel('Time'); ylabel('Stretched Depth [z_m = int Ndz / N_m]');
c = colorbar; c.Label.String = 'Shear [1/s]';
title('dv/dz')




%f = figure(); clf; set(gcf,'Color','w'); hold on 
ax(1) = subplot(2,2,1); hold on 
p1 = pcolor(dates(idx0:idx1),z_uv,shear_u(:,idx0:idx1)); %,'DisplayName','Zonal Velocity'); %'Color',orange);
xlabel('Time'); ylabel('Depth [m]');
shading flat; grid on; set(gca, 'GridLineStyle', '--','GridAlpha', 0.7,'Layer','top') 
set(gca, 'YDir', 'reverse');
c = colorbar; c.Label.String = 'Shear [1/s]';
title('du/dz')


ax(2) = subplot(2,2,2); hold on 
p2 = pcolor(dates(idx0:idx1),z_uv,shear_v(:,idx0:idx1)); %,'DisplayName','Meridional Velocity'); %,'Color',orange);
shading flat; grid on; set(gca, 'GridLineStyle', '--','GridAlpha', 0.7,'Layer','top'); % Slight transparency
set(gca, 'YDir', 'reverse');
xlabel('Time'); ylabel('Depth [m]');
c = colorbar; c.Label.String = 'Shear [1/s]';
title('dv/dz')
