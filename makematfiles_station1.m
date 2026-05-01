% script for concatenating and processing u/v velocities from SLE stations 

% data_dir = /change this to kerstin's harddrive ; 
data_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/processed_ship/'; 
out_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/processed_ship/';

station = 1; 
ladcp_dir = [data_dir 'ladcp/station' num2str(station) '/']; 

cast = 1; 
file = sprintf('%02d.mat',cast);
fn = [ladcp_dir file]; 
S = load(fn); 

% load variables
u_cast = S.dr.u; 
v_cast = S.dr.v; 

lon_cast = S.dr.lon; 
lat_cast = S.dr.lat; 

p_cast = S.dr.p; 

% check the uv ctd 
u_ctd = S.dr.uctd; 
v_ctd = S.dr.vctd; 
w_ctd = S.dr.wctd; 

% for TS
t_ctd = S.dr.ctd_t;
s_ctd = S.dr.ctd_s;

% separate out up cast and downcast 

u_do = S.dr.u_do;
v_do = S.dr.v_do; 

u_up = S.dr.u_up; 
v_up = S.dr.v_up; 


%% process and concatenate 

% for station 1: 18 casts 
if station == 1 
    casts = 1:18; 
elseif station == 2
    casts = 1:36; 
end
    
numcasts = length(casts);

% preallocate for storage

if station == 1
    zlevels = 559; 
    ctdlevs = 838; 
elseif station == 2
    zlevels = 234; 
    ctdlevs = 328; 
end


u_all = nan(numcasts,zlevels); 
v_all = nan(numcasts,zlevels); 
p_all = nan(numcasts,zlevels); 

u_down = nan(numcasts,zlevels); 
v_down = nan(numcasts,zlevels); 
p_down = nan(numcasts,zlevels); 
time_down = NaT(numcasts,1);


% separate out upcast and downcast
u_grid = nan(numcasts*2,zlevels); 
v_grid = nan(numcasts*2,zlevels); 
p_grid = nan(numcasts*2,zlevels); 
time_grid = NaT(numcasts*2,zlevels); 
time_along = NaT(numcasts,ctdlevs); 
time_updo = NaT(numcasts*2,1);

u_shearmethod = nan(numcasts,zlevels); 
v_shearmethod = nan(numcasts,zlevels); 
w_shearmethod = nan(numcasts,zlevels); 


t_all = nan(numcasts,zlevels); 
s_all = nan(numcasts,zlevels); 

u_ctd = nan(numcasts,ctdlevs); 
v_ctd = nan(numcasts,ctdlevs); 
w_ctd = nan(numcasts,ctdlevs); 

% ship info 
ship_lon = nan(numcasts,ctdlevs); 
ship_lat = nan(numcasts,ctdlevs); 

ladcp_lon = nan(numcasts,1); 
ladcp_lat = nan(numcasts,1); 
ladcp_time = NaT(numcasts,1);


zlevcount = []; 
ctdlevcount = []; 

for cc = casts
 
    if cc == 21
        % skip cast 21 still  
        % leave as nan 

    else    
        file = sprintf('%02d.mat',cc);
        fn = [ladcp_dir file]; 
        S = load(fn); 
        
        % load variables
        u_cast = S.dr.u; 
        v_cast = S.dr.v; 
    
        u_goingdo = S.dr.u_do; v_goingdo = S.dr.v_do; 
        u_goingup = S.dr.u_up; v_goingup = S.dr.v_up; 
    
        time_start = datetime(S.dr.tim(1),'ConvertFrom','juliandate'); 
        time_end = datetime(S.dr.tim(end),'ConvertFrom','juliandate'); 
        time_middle = time_start + (time_end-time_start)/2 ; 
    
        ushear = S.dr.u_shear_method; 
        vshear = S.dr.v_shear_method; 
        wshear = S.dr.w_shear_method; 
        
        
        lon_cast = S.dr.lon; 
        lat_cast = S.dr.lat; 
        
        p_cast = S.dr.p; 
        date_cast = datetime(S.dr.date); 
    
        disp(length(p_cast))
        zlevcount = cat(1,zlevcount,length(p_cast));
        
        % check the uv ctd 
        u_ctd_cast = S.dr.uctd; 
        v_ctd_cast = S.dr.vctd; 
        w_ctd_cast = S.dr.wctd; 
        
        % ship info 
    
        slon = S.dr.shiplon; 
        slat = S.dr.shiplat; 
    
        disp(length(u_ctd_cast))
        ctdlevcount = cat(1,ctdlevcount,length(u_ctd_cast)) ;
    
        % save T, S - though bad data 
        t_ctd = S.dr.ctd_t;
        s_ctd = S.dr.ctd_s;
    
        % concatenate 
    
        ladcp_lon(cc) = lon_cast; 
        ladcp_lat(cc) = lat_cast;
        ladcp_time(cc) = date_cast;
    
        u_all(cc,1:length(p_cast)) = u_cast; 
        v_all(cc,1:length(p_cast)) = v_cast; 
        p_all(cc,1:length(p_cast)) = p_cast; 

        u_down(cc,1:length(p_cast)) = u_goingdo; 
        v_down(cc,1:length(p_cast)) = v_goingdo; 
        p_down(cc,1:length(p_cast)) = p_cast; 
        time_down(cc) = time_start;

    
        % downcasts
        u_grid(1+(cc-1)*2,1:length(p_cast)) = u_goingdo; 
        v_grid(1+(cc-1)*2,1:length(p_cast)) = v_goingdo; 
        p_grid(1+(cc-1)*2,1:length(p_cast)) = p_cast; 
    
        % downcast aligned with cast start time 
        time_grid(1+(cc-1)*2,1:length(p_cast)) = time_start; 
        time_updo(1+(cc-1)*2) = time_start; 
        
    
        % upcasts
        u_grid(cc*2,1:length(p_cast)) = u_goingup; 
        v_grid(cc*2,1:length(p_cast)) = v_goingup; 
        p_grid(cc*2,1:length(p_cast)) = p_cast; 
    
        % upcast aligned with cast mid time 
        time_grid(cc*2,1:length(p_cast)) = time_middle; 
        time_updo(cc*2) = time_middle; 
    
    
        u_shearmethod(cc,1:length(p_cast)) = ushear; 
        v_shearmethod(cc,1:length(p_cast)) = vshear; 
        w_shearmethod(cc,1:length(p_cast)) = wshear; 
    
    
        t_all(cc,1:length(p_cast)) = t_ctd; 
        s_all(cc,1:length(p_cast)) = s_ctd; 
        
        u_ctd(cc,1:length(u_ctd_cast)) = u_ctd_cast; 
        v_ctd(cc,1:length(u_ctd_cast)) = v_ctd_cast;
        w_ctd(cc,1:length(u_ctd_cast)) = w_ctd_cast;
    
        ship_lon(cc,1:length(u_ctd_cast)) = slon;    
        ship_lat(cc,1:length(u_ctd_cast)) = slat;
    
        time_along(cc,1:length(u_ctd_cast)) = datetime(S.dr.tim,'ConvertFrom','juliandate');

    end

end

%% now plot 
addpath /home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/st_lucia_analysis

cast_dim = repmat(casts,zlevels,1)'; 
time_dim = repmat(ladcp_time,1,zlevels); 

figure()
subplot(3,2,1); hold on 
title('u velocity')
pcolor(time_dim,p_all,u_all); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
xtickformat('dd HH:mm')
set(gca, 'YDir', 'reverse')
cb = colorbar; cmocean('balance')
ylabel(cb,'u [m/s]','Rotation',270)
clim([-0.2 0.2]); 
%xlim([datetime(2025,2,28,13,0,0) datetime(2025,3,1,15,0,0)])

subplot(3,2,2); hold on 
title('v velocity')
pcolor(time_dim,p_all,v_all); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
cb = colorbar; cmocean('balance')
ylabel(cb,'v [m/s]','Rotation',270)
clim([-0.2 0.2]); 
%xlim([datetime(2025,2,28,13,0,0) datetime(2025,3,1,15,0,0)])

subplot(3,2,3); hold on 
title('u velocity: up/down')
pcolor(time_grid,p_grid,u_grid); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
xtickformat('dd HH:mm')
set(gca, 'YDir', 'reverse')
cb = colorbar; cmocean('balance')
ylabel(cb,'u [m/s]','Rotation',270)
clim([-0.2 0.2]); 
%xlim([datetime(2025,2,28,13,0,0) datetime(2025,3,1,15,0,0)])

subplot(3,2,4); hold on 
title('v velocity: up/down')
pcolor(time_grid,p_grid,v_grid); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
cb = colorbar; cmocean('balance')
ylabel(cb,'v [m/s]','Rotation',270)
clim([-0.2 0.2]); 
%xlim([datetime(2025,2,28,13,0,0) datetime(2025,3,1,15,0,0)])


subplot(3,2,5); hold on 
title('u velocity: down only')
pcolor(time_dim,p_down,u_down); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
xtickformat('dd HH:mm')
set(gca, 'YDir', 'reverse')
cb = colorbar; cmocean('balance')
ylabel(cb,'u [m/s]','Rotation',270)
clim([-0.2 0.2]); 
%xlim([datetime(2025,2,28,13,0,0) datetime(2025,3,1,15,0,0)])

subplot(3,2,6); hold on 
title('v velocity: down only')
pcolor(time_dim,p_down,v_down); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
cb = colorbar; cmocean('balance')
ylabel(cb,'v [m/s]','Rotation',270)
clim([-0.2 0.2]); 
%xlim([datetime(2025,2,28,13,0,0) datetime(2025,3,1,15,0,0)])


%% calculate shear 


[nc,nz] = size(p_grid);
shear_u = nan(nc,nz);
shear_v = nan(nc,nz);

disp(size(shear_u))
for cc = 1:nc 
    % need to do individually for each column because they each start at a different point due to NaNs
    
    for ii = 1:nz-1
        % going from the top to bottom  
        
        if isnan(u_grid(cc,ii+1)) == 0
            % do shear calculations only if the level below has a value 
            shear_u(cc,ii) = (u_grid(cc,ii+1)-u_grid(cc,ii)) / (p_grid(cc,ii+1) - p_grid(cc,ii)); 
            shear_v(cc,ii) = (v_grid(cc,ii+1)-v_grid(cc,ii)) / (p_grid(cc,ii+1) - p_grid(cc,ii));

        end 
    end
end


%% now plot shear

figure()
subplot(2,1,1); hold on 
title('u shear')
pcolor(time_grid,p_grid,shear_u); 
shading flat 
xlabel('Cast #'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cb = colorbar; cmocean('balance')
ylabel(cb,'du/dz [1/s]','Rotation',270)
clim([-0.015 0.015])

subplot(2,1,2); hold on 
title('v shear')
pcolor(time_grid,p_grid,shear_v); 
shading flat 
xlabel('Cast #'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cb = colorbar; cmocean('balance')
ylabel(cb,'dv/dz [1/s]','Rotation',270)
clim([-0.015 0.015])

%% based on wavenumber analysis, need to interpolate to 20m to calculate shear

% can we do a rolling smoother or do we need to straight up coarsen? 
% lets try both 


[nc,nz] = size(u_grid);
rolling_shear_u = nan(nc,nz);
rolling_shear_v = nan(nc,nz);
nz_20m = 5; % roughly 5 datapoints in 20m 

bottom_depth = max(p_grid,[],'all'); 
interp_grid = 0:20:bottom_depth; 
interp_nz = length(interp_grid);
shear_grid = nan(1,interp_nz-1);

interp_u = nan(nc,interp_nz);
interp_v = nan(nc,interp_nz);

interp_shear_u = nan(nc,interp_nz-1);
interp_shear_v = nan(nc,interp_nz-1);
tot_shear = nan(nc,interp_nz-1);


for cc = 1:nc 
    % need to do individually for each column because they each start at a different point due to NaNs
    
    % rolling method 
    deepest = sum(~isnan(p_grid(cc,:)));
    for ii = 3:deepest-3

        % going from the top to bottom - calculate rolling
        % 5 datapoints in 20m --> take 2 on each side 
        
        u_i = mean(u_grid(cc,ii-2:ii+2),'omitnan'); 
        v_i = mean(v_grid(cc,ii-2:ii+2),'omitnan'); 
        p_i = mean(p_grid(cc,ii-2:ii+2),'omitnan');
        u_ii = mean(u_grid(cc,ii-1:ii+3),'omitnan');
        v_ii = mean(v_grid(cc,ii-1:ii+3),'omitnan');
        p_ii = mean(p_grid(cc,ii-1:ii+3),'omitnan');

        rolling_shear_u(cc,ii) = ( u_ii-u_i) / (p_ii-p_i ); 
        rolling_shear_v(cc,ii) = ( v_ii-v_i) / (p_ii-p_i ); 
        
    end

    % coarsening method 
    for ii = 1:interp_nz
        depth_on_grid = interp_grid(ii);
        [~,depth_index] = min(abs(p_grid(1,:)-depth_on_grid));

        if ii == 1
            % we are at the surface - only use top 2 in the mean 
            interp_u(cc,ii) = mean(u_grid(cc,depth_index:depth_index+2),'omitnan');
            interp_v(cc,ii) = mean(v_grid(cc,depth_index:depth_index+2),'omitnan');

        elseif ii == interp_nz
            % we are at the bottom - only use the bottom available in the mean
            interp_u(cc,ii) = mean(u_grid(cc,depth_index-2:end),'omitnan');
            interp_v(cc,ii) = mean(v_grid(cc,depth_index-2:end),'omitnan');

        else 
            interp_u(cc,ii) = mean(u_grid(cc,depth_index-2:depth_index+2),'omitnan');
            interp_v(cc,ii) = mean(v_grid(cc,depth_index-2:depth_index+2),'omitnan');

        end

    end 

    % now, go through that loop but to calculate shear 
    for ii = 1:interp_nz-1
        shear_grid(ii) = (interp_grid(ii+1) + interp_grid(ii)) / 2 ;
        interp_shear_u(cc,ii) = (interp_u(cc,ii+1)-interp_u(cc,ii)) / (interp_grid(ii+1) - interp_grid(ii)); 
        interp_shear_v(cc,ii) = (interp_v(cc,ii+1)-interp_v(cc,ii)) / (interp_grid(ii+1) - interp_grid(ii)); 

        tot_shear(cc,ii) = interp_shear_u(cc,ii)^2 + interp_shear_v(cc,ii)^2;
    end
end

[x,y] = meshgrid(time_updo,interp_grid);
[xshear,yshear] = meshgrid(time_updo,shear_grid); 

fig1=figure()
subplot(3,2,1); hold on 
[t,s] = title('Station 1','interpolated u velocity');
pcolor(x',y',interp_u); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'u velocity [m/s]','Rotation',270)
clim([-0.25 0.25])

subplot(3,2,2); hold on 
[t,s] = title('Station 1','interpolated v velocity');
pcolor(x',y',interp_v); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'v velocity [m/s]','Rotation',270)
clim([-0.25 0.25])


subplot(3,2,3); hold on 
subtitle('interpolated u shear')
pcolor(xshear',yshear',interp_shear_u); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'du/dz [1/s]','Rotation',270)
clim([-0.01 0.01])


subplot(3,2,4); hold on 
subtitle('interpolated v shear')
pcolor(xshear',yshear',interp_shear_v); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'dv/dz [1/s]','Rotation',270)
clim([-0.01 0.01])

subplot(3,2,5); hold on 
subtitle('total shear')
pcolor(xshear',yshear',tot_shear); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
set(gca,'ColorScale','log');
ylabel(cb,'(du/dz )^2 + (dv/dz)^2 [1/s2]','Rotation',270)

%%
subplot(3,2,5); hold on 
subtitle('rolling u shear')
pcolor(time_grid, p_grid,rolling_shear_u); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'du/dz [1/s]','Rotation',270)
clim([-0.01 0.01])

subplot(3,2,6); hold on 
subtitle('rolling v shear')
pcolor(time_grid, p_grid,rolling_shear_v); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'dv/dz [1/s]','Rotation',270)
clim([-0.01 0.01])



%% first pass variable save 

% station 1 

LADCP1 = struct; 
LADCP1.lon = lon; LADCP1.lat = lat; LADCP1.time = time; 
LADCP1.u = u_all; LADCP1.v = v_all; LADCP1.p = p_all; 
LADCP1.u_shearmethod = u_shearmethod; LADCP1.v_shearmethod = v_shearmethod; LADCP1.w_shearmethod = w_shearmethod; 
LADCP1.tctd = t_ctd; LADCP1.sctd = s_ctd; 
LADCP1.u_ctd = u_ctd; LADCP1.v_ctd = v_ctd; LADCP1.w_ctd = w_ctd; 
LADCP1.ship_lon = ship_lon; LADCP1.ship_lat = ship_lat; 

save([out_dir 'LADCP_station' num2str(station)],'LADCP1'); 

%%
% save ladcp depths to interpolate with CTD 
ladcp_p = p_all; 

%% station 1 struct 

% what to save: 
% time dimension (at each grid point); time coordinate (for each cast)
% depth dimension; depth coordinate

% raw u/v separated into up/down. gridded (20m) u/v. only down u/v
% gridded shear (20m) 

STAT1 = struct; 
STAT1.ladcp_lon = ladcp_lon; STAT1.ladcp_lat = ladcp_lat; 
STAT1.ladcp_time_coord = time_updo; STAT1.ladcp_time_grid = time_grid; 
STAT1.u = u_grid; STAT1.v = v_grid; 
STAT1.ladcp_p_coord = interp_grid; STAT1.ladcp_p_grid = p_grid; 

STAT1.interp_u = interp_u; STAT1.interp_v = interp_v; 

STAT1.shear_coord = shear_grid; 
STAT1.interp_shear_u = interp_shear_u; STAT1.interp_shear_v = interp_shear_v; 


%% now use kerstin's code to correct temp and salinity from t1/s1  vs. t2/s2 

% based on those plots this is a first best guess on which is "prettiest"
goodT = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1];
goodS = [1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1];
 

casts = 1:18; 

% set up matrices for putting the data in
Temp1 = nan(max(lengths),18);
Sal1 = Temp1; Pres1 = Temp1;
dn1 = nan(1,18);

for ii = casts
    load(filelist{ii});

    if goodT(ii) == 1
    T = datad.t1;
    else
    T = datad.t2;
    end

    if goodS(ii) == 1
    S = datad.s1;
    else
    S = datad.s2;
    end

    Temp1(1:lengths(ii),ii) = fillmissing(T,'linear');
    Sal1(1:lengths(ii),ii) = fillmissing(S,'linear');
    Pres1(1:lengths(ii),ii) = fillmissing(datad.p,'linear');
    dn1(ii) = nanmean(datad.datenum);
   
end


%% next, concatenate temp and salinity from .mat files 

ctd_dir = [data_dir 'ctd_mat/station' num2str(station) '/']; 

station = 1;
cruise = 'SR2503';

if station == 1
    ctd_dir = [data_dir 'ctd_mat/station' num2str(station) '/' ];
elseif station == 2
    ctd_dir = [data_dir 'ctd_mat/station' num2str(station) '/station_2' ];
end


cast = 1; 
file = sprintf('%s_POstation_%d_cast_%02d.mat',cruise,station,cast);
fn = [ctd_dir file]; 
S = load(fn); 

t1d = S.datad.t1; 
s1d = S.datad.s1; 
t2d = S.datad.t2; 
s2d = S.datad.s2; 

t1u = S.datau.t1; 
s1u = S.datau.s1; 
t2u = S.datau.t2; 
s2u = S.datau.s2; 


theta1 = S.datad.theta1; 
sigma1 = S.datad.sigma1; 
theta2 = S.datad.theta2; 
sigma2 = S.datad.sigma2; 

lon = S.datad.lon; 
lat = S.datad.lat; 

depth = S.datad.depth; 
p = S.datad.p; 
date = datetime(S.datad.datenum,'ConvertFrom','datenum'); 
time = S.datad.time; 

%% initial plot 
datadir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/processed_ship/ctd_mat/station1/';
files = dir(fullfile(datadir, 'SR2503_POstation_1_cast_*.mat'));
filelist = fullfile({files.folder}, {files.name})';


% first a visual inspection
for ii = 1:18
    load(filelist{ii});

    t1 = datad.t1;
    t2 = datad.t2;
    s1 = datad.s1;
    s2 = datad.s2;

    lengths(ii) = length(t1);

    figure(ii)
    clf
    subplot(1,2,1)
    plot([t1 t2],datad.depth)
    legend('t1','t2'); axis ij
    xlim([2 15])
    subplot(1,2,2)
    plot([s1 s2],datad.depth)
    legend('s1','s2'); axis ij
    xlim([30 36])
    title(['Cast number: ' num2str(ii)])
    
end
clear datad datau datad_1m datau_1m

%% process and concatenate 

addpath /home/vboatwright/OneDrive/Documents/SIO/projects/gsw_matlab
addpath /home/vboatwright/OneDrive/Documents/SIO/projects/gsw_matlab/library
addpath /home/vboatwright/OneDrive/Documents/SIO/projects/gsw_matlab/thermodynamics_from_t
addpath /home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/st_lucia_analysis


% for station 1: 18 casts 
if station == 1 
    casts = 1:18; 
elseif station == 2
    casts = 1:36; 
end
    
numcasts = length(casts);

% preallocate based on cast 1 

array_sizes = []; 

if station == 1
    zlevels = 8921; 
elseif station == 2
    zlevels = 3729; 
end

% separate out upcast and downcast - will decide per cast which sensor is best
temp = nan(numcasts*2,zlevels); 
sal = nan(numcasts*2,zlevels); 
p = nan(numcasts*2,zlevels); 
rho = nan(numcasts*2,zlevels); 
n2 = nan(numcasts*2,zlevels); 
timectd = nan(numcasts*2,zlevels);

temp_down = nan(numcasts,zlevels); 
sal_down = nan(numcasts,zlevels); 
p_down = nan(numcasts,zlevels); 
rho_down = nan(numcasts,zlevels); 
n2_down = nan(numcasts,zlevels); 
time_down = nan(numcasts,zlevels); 

t_grid = nan(numcasts*2,zlevels); 
s_grid = nan(numcasts*2,zlevels); 
p_grid = nan(numcasts*2,zlevels); 
date_grid = NaT(numcasts*2,zlevels); 
time_grid = nan(numcasts*2,zlevels); 

lon_grid = nan(numcasts*2,zlevels); 
lat_grid = nan(numcasts*2,zlevels); 

time_coord = nan(numcasts*2,1);

goodT = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,1,1,1];
goodS = [1,1,1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,1];

for cc = casts
   
    disp(cc)
    if cc == 21 
        % no file  
        % leave as nan 

    else

    file = sprintf('%s_POstation_%d_cast_%02d.mat',cruise,station,cc);
    disp(file)
    fn = [ctd_dir file]; 
    S = load(fn); 

    % load variables
    zlength = length(S.datad.t1); 

    % temp 

    if goodT(ii) == 1
        % down going 
        t_goingdo = S.datad.t1; 
        % up going 
        t_goingup = S.datau.t1; 

    else
        % down going 
        t_goingdo = S.datad.t2; 
        % up going 
        t_goingup = S.datau.t2; 
    end

    % salinity 

    if goodS(ii) == 1
        % down going 
        s_goingdo = S.datad.s1; 
        % up going 
        s_goingup = S.datau.s1; 

    else
        % down going 
        s_goingdo = S.datad.s2; 
        % up going 
        s_goingup = S.datau.s2; 
    end


    if cc == 10
        sfilled = fillmissing(s_goingdo,'linear'); 
        s_goingdo = sfilled; 
        sfilled = fillmissing(s_goingup,'linear'); 
        s_goingup = sfilled; 
    end

    if cc == 12
        sfilled = fillmissing(s_goingdo,'linear'); 
        s_goingdo = sfilled; 
        sfilled = fillmissing(s_goingup,'linear'); 
        s_goingup = sfilled;     
    end


    % other vars 
    % down going 
    depth_down = S.datad.depth; 
    p_goingdo = S.datad.p; 
    date_goingdo = datetime(S.datad.datenum,'ConvertFrom','datenum'); 
    time_goingdo = S.datad.time; 

    lon_goingdo = S.datad.lon; 
    lat_goingdo = S.datad.lat; 
    
    array_sizes = cat(1,array_sizes,zlength); 

    % other vars
    % up going 
    depth_up = S.datau.depth; 
    p_goingup = S.datau.p; 
    date_goingup = datetime(S.datau.datenum,'ConvertFrom','datenum'); 
    time_goingup = S.datau.time; 
    lon_goingup = S.datau.lon; 
    lat_goingup = S.datau.lat; 


    % GSW to down 
    sa_do = gsw_SA_from_SP(s_goingdo,p_goingdo,239,34);
    ct_do = gsw_CT_from_t(sa_do,t_goingdo,p_goingdo);
    rho_do = gsw_rho(sa_do,ct_do,p_goingdo)-1000;
    rho_do_sort = sort(rho_do);
    N2_do_sort = smoothdata(smoothdata(9.8./nanmean(rho_do_sort).*diffdiff(sort(rho_do_sort),1)*2,'gaussian',32),2,'gaussian',5);

    % GSW to up 
    sa_up = gsw_SA_from_SP(s_goingup,p_goingup,239,34);
    ct_up = gsw_CT_from_t(sa_up,t_goingup,p_goingup);
    rho_up = gsw_rho(sa_up,ct_up,p_goingup)-1000;
    rho_up_sort = sort(rho_up);
    N2_up_sort = smoothdata(smoothdata(9.8./nanmean(rho_up_sort).*diffdiff(sort(rho_up_sort),1)*2,'gaussian',32),2,'gaussian',5);


    % concatenate 

    temp_down(cc,1:zlength) = ct_do; 
    sal_down(cc,1:zlength) = sa_do; 
    p_down(cc,1:zlength) = p_goingdo; 
    rho_down(cc,1:zlength) = rho_do_sort; 
    n2_down(cc,1:length(N2_do_sort)) = N2_do_sort;  
    time_down(cc) = mean(time_goingdo(find(~isnan(time_goingdo))));

    % put it all together on the _grid files: 


    % downcasts
    temp(1+(cc-1)*2,1:zlength) = ct_do; 
    sal(1+(cc-1)*2,1:zlength) = sa_do; 
    p(1+(cc-1)*2,1:zlength) = p_goingdo;  
    rho(1+(cc-1)*2,1:zlength) = rho_do_sort;  
    n2(1+(cc-1)*2,1:length(N2_do_sort)) = N2_do_sort;  
    date_grid(1+(cc-1)*2,1:zlength) = date_goingdo;     
    time_grid(1+(cc-1)*2,1:zlength) = time_goingdo;     

    lon_grid(1+(cc-1)*2,1:zlength) = lon_goingdo; 
    lat_grid(1+(cc-1)*2,1:zlength) = lat_goingdo; 

    time_coord(1+(cc-1)*2) = mean([S.datad.time,S.datad.time],'all','omitnan');

    % upcasts
    temp(cc*2,1:zlength) = ct_up; 
    sal(cc*2,1:zlength) = sa_up; 
    p(cc*2,1:zlength) = p_goingup; 
    rho(cc*2,1:zlength) = rho_up_sort; 
    n2(cc*2,1:length(N2_up_sort)) = N2_up_sort; 
    date_grid(cc*2,1:zlength) = date_goingup;
    time_grid(cc*2,1:zlength) = time_goingup;     

    lon_grid(cc*2,1:zlength) = lon_goingup; 
    lat_grid(cc*2,1:zlength) = lat_goingup; 

    time_coord(cc*2) = mean([S.datau.time,S.datau.time],'all','omitnan');

    % gridded?
    % time_updo(cc*2) = time_middle; 
    end

end


%% now plot 

figure()
subplot(3,1,1); hold on 
title('cons temp')
pcolor(date_grid,p,temp); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
xtickformat('dd HH:mm')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('thermal')); 
cb = colorbar;
ylabel(cb,'[deg C]','Rotation',270)
% clim([-0.2 0.2])

subplot(3,1,2); hold on 
title('sal')
pcolor(date_grid,p,sal); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('haline')); 
cb = colorbar;
ylabel(cb,'[psu]','Rotation',270)
%clim([32 36])

subplot(3,1,3); hold on 
title('rho')
pcolor(date_grid,p,rho); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('dense')); 
cb = colorbar;
ylabel(cb,'[psu]','Rotation',270)
%clim([32 36])


%% but for now just save designated values 

% station 1 

CTD1 = struct; 
CTD1.lons = lons; CTD1.lats = lats; CTD1.times = times; CTD1.dates = dates; 
CTD1.t1 = t1_all; CTD1.s1 = s1_all; CTD1.p = p_all; 
CTD1.t2 = t2_all; CTD1.s2 = s2_all; CTD1.depth = depth_all; 

save([out_dir 'CTD_station' num2str(station)],'CTD1'); 

%% interpolate to coarse grid to save 

% interpolate CTD to u/v grid 
% use: interp_grid and nc, interp_nz

interp_t = nan(nc,interp_nz);
interp_s = nan(nc,interp_nz);
interp_rho = nan(nc,interp_nz); 

interp_n2 = nan(nc,interp_nz-1);

for cc = 1:nc 
    % need to do individually for each column because they each start at a different point due to NaNs
    disp(cc)
    % coarsening method 
    for ii = 1:interp_nz
        depth_on_grid = interp_grid(ii);
        [~,depth_index] = min(abs(p(1,:)-depth_on_grid));
        [~,shallowidx] = min(abs(p(1,:)-(depth_on_grid-10)));
        [~,deepidx] = min(abs(p(1,:)-(depth_on_grid+10)));

        % we are at the surface - only use top 2 in the mean 
        interp_t(cc,ii) = mean(temp(cc,shallowidx:deepidx),'omitnan');
        interp_s(cc,ii) = mean(sal(cc,shallowidx:deepidx),'omitnan');
        interp_rho(cc,ii) = mean(rho(cc,shallowidx:deepidx),'omitnan');

    end

    [cast_n2, cast_n2_grid] = gsw_Nsquared(interp_s(cc,:),interp_t(cc,:),interp_grid,mean(ladcp_lat));
    interp_n2(cc,:) = cast_n2;
    n2_grid = cast_n2_grid; % each of these should be the same 

end 

% now calculate n2 
% [interp_N2, n2_grid] = gsw_Nsquared(interp_s',interp_t',interp_grid',mean(ladcp_lat));

%%


[X,Y] = meshgrid(time_coord,interp_grid);
[XN2,YN2] = meshgrid(time_coord,n2_grid); 

fig1=figure()
subplot(2,2,1); hold on 
[t,s] = title('Station 1','interpolated conservative temp');
pcolor(X',Y',interp_t); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('thermal'); cb = colorbar;
ylabel(cb,'temp [degC]','Rotation',270)

subplot(2,2,2); hold on 
[t,s] = title('Station 1','interpolated absolute salinity');
pcolor(X',Y',interp_s); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('haline'); cb = colorbar;
ylabel(cb,'salinity','Rotation',270)


subplot(2,2,3); hold on 
subtitle('interpolated density')
pcolor(X',Y',interp_rho); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('dense'); cb = colorbar;
ylabel(cb,'rho [kg/m3]','Rotation',270)


subplot(2,2,4); hold on 
subtitle('interpolated n2')
pcolor(XN2',YN2',interp_n2); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('matter'); cb = colorbar;
set(gca,'ColorScale','log');
ylabel(cb,'N2 [1/s2]','Rotation',270)
%clim([-0.01 0.01])








%% add CTD to STAT1 struct 

STAT1.ctd_time_coord = time_coord; STAT1.ctd_date_grid = date_grid; 
STAT1.t = temp; STAT1.s = sal; 
STAT1.ctd_p = p; 

STAT1.ctd_time_do = time_down; STAT1.n2_down = n2_down; 
STAT1.t_do = temp_down; STAT1.s_do = sal_down; 
STAT1.ctd_p_do = p_down; 
  
STAT1.interp_t = interp_t; STAT1.interp_s = interp_s; 

STAT1.interp_rho = interp_rho; 
STAT1.n2_coord = n2_grid; 
STAT1.interp_n2 = interp_n2;


%% calculate richardson 

% Ri = N2 / S2 
% mixing at high ri-1 
% N2 = -g/rho d rho / dz 
% [N2, p_mid] = gsw_Nsquared(SA,CT,p,{lat})

inv_Ri = tot_shear./interp_n2;  
Ri = interp_n2./tot_shear;  


%% now plot it! 

figure(); hold on 
subplot(3,1,1); hold on 
[t,s] = title('Station 1','Interpolated N2');
pcolor(xshear',yshear',interp_n2); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('matter'); cb = colorbar;
set(gca,'ColorScale','log');
ylabel(cb,'N2 [1/s2]','Rotation',270)


subplot(3,1,2); hold on 
title('Total Shear Squared');
pcolor(xshear',yshear',tot_shear); 
shading interp 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('matter'); cb = colorbar;
set(gca,'ColorScale','log');
ylabel(cb,'Shear [1/s2]','Rotation',270)


peaks = [XN2',YN2',inv_Ri]; 

subplot(3,1,3); hold on
subtitle('Inverse Richardson')
pcolor(xshear',yshear',inv_Ri); 
shading interp
% contour(peaks,[0.001,1],'Color','white'); 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('thermal'); cb = colorbar;
set(gca,'ColorScale','log');
ylabel(cb,'1/Ri','Rotation',270)


%%
% 
xticks(time_on_grid(1:2:end,1))
xticklabels(datestr(time_on_grid(1:2:end,1),'mm/dd HH:MM'))


timenums = datenum(time_on_grid);
clevels = [1];

subplot(2,1,2); hold on
pcolor(timenums(:,1:end-1),p_gridded(:,1:end-1),inv_Ri')
contour(timenums(:,1:end-1),p_gridded(:,1:end-1),inv_Ri',clevels,'Color','k')
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('matter')); 
cb = colorbar;
ylabel(cb,'Ri^{-1} #','Rotation',270)
set(gca,'ColorScale','log')
xticks(timenums(1:2:end,1))
xticklabels(datestr(time_on_grid(1:2:end,1),'mm/dd HH:MM'))
























%% now process the 1m binned data
% coming soon:  of both down and up data!!! 

cast = 1; 
file = sprintf('%s_POstation_%d_cast_%02d.mat',cruise,station,cast);
fn = [ctd_dir file]; 
S = load(fn); 

t1 = S.datad_1m.t1; 
s1 = S.datad_1m.s1; 
t2 = S.datad_1m.t2; 
s2 = S.datad_1m.s2; 


theta1 = S.datad_1m.theta1; 
sigma1 = S.datad_1m.sigma1; 
theta2 = S.datad_1m.theta2; 
sigma2 = S.datad_1m.sigma2; 

lon = S.datad_1m.lon; 
lat = S.datad_1m.lat; 

depth = S.datad_1m.depth; 
p = S.datad_1m.p; 
date = datetime(S.datad_1m.datenum,'ConvertFrom','datenum'); 
time = S.datad_1m.time; 

%% initial plot 

figure()
plot(t1,depth); hold on 
plot(s1,depth)

%% process 1m data and concatenate 

% for station 1: 18 casts 
if station == 1 
    casts = 1:18; 
elseif station == 2
    casts = 1:36; 
end
    
numcasts = length(casts);

% preallocate based on cast 1's 1m bins 

array_sizes = []; 

if station == 1
    zlevels = 2231; 
elseif station == 2
    zlevels = 0; % havent set this yet  
end

t1_all = nan(numcasts,zlevels); 
s1_all = nan(numcasts,zlevels); 
p_all = nan(numcasts,zlevels); 

t2_all = nan(numcasts,zlevels); 
s2_all = nan(numcasts,zlevels); 
depth_all = nan(numcasts,zlevels); 


% separate out upcast and downcast
t1_grid = nan(numcasts*2,zlevels); 
s1_grid = nan(numcasts*2,zlevels); 
cp_grid = nan(numcasts*2,zlevels); 
ctime_grid = NaT(numcasts*2,zlevels); 

ctd_lons = nan(numcasts,zlevels); 
ctd_lats = nan(numcasts,zlevels); 
ctd_dates = NaT(numcasts,zlevels); 
ctd_times = nan(numcasts,zlevels); 

for cc = casts
   
    disp(cc)
    if cc == 21 
        % no file  
        % leave as nan 

    else

    file = sprintf('%s_POstation_%d_cast_%02d.mat',cruise,station,cc);
    disp(file)
    fn = [ctd_dir file]; 
    S = load(fn); 

    
    % load variables
    
    t1cast = S.datad_1m.t1; 
    s1cast = S.datad_1m.s1; 
    t2cast = S.datad_1m.t2; 
    s2cast = S.datad_1m.s2; 


    t1down = S.datad_1m.t1; 
    s1down = S.datad_1m.s1; 
    t2down = S.datad_1m.t2; 
    s2down = S.datad_1m.s2; 
    pdown = S.datad_1m.p; 

    t1up = S.datau_1m.t1; 
    s1up = S.datau_1m.s1; 
    t2up = S.datau_1m.t2; 
    s2up = S.datau_1m.s2; 
    pup = S.datau_1m.p; 

    theta1cast = S.datad_1m.theta1; 
    sigma1cast = S.datad_1m.sigma1; 
    theta2cast = S.datad_1m.theta2; 
    sigma2cast = S.datad_1m.sigma2; 
    
    loncast = S.datad_1m.lon; 
    latcast = S.datad_1m.lat; 
    
    depthcast = S.datad_1m.depth; 
    pcast = S.datad_1m.p; 
    disp(length(pcast))
    datecast = datetime(S.datad_1m.datenum,'ConvertFrom','datenum'); 
    timecast = S.datad_1m.time; 

    array_sizes = cat(1,array_sizes,length(pcast)); 

    % concatenate 

    t1_all(cc,1:length(pcast)) = t1cast; 
    s1_all(cc,1:length(pcast)) = s1cast; 
    p_all(cc,1:length(pcast)) = pcast; 

    t2_all(cc,1:length(pcast)) = t2cast; 
    s2_all(cc,1:length(pcast)) = s2cast; 
    depth_all(cc,1:length(pcast)) = depthcast; 

    time_down = datetime(S.datad_1m.time,'ConvertFrom','posixtime'); 
    time_up = datetime(S.datau_1m.time,'ConvertFrom','posixtime'); 

    % concatenate with separated up and down 

    % downcasts
    t1_grid(1+(cc-1)*2,1:length(pcast)) = t1down; 
    s1_grid(1+(cc-1)*2,1:length(pcast)) = s1down; 
    cp_grid(1+(cc-1)*2,1:length(pcast)) = pdown; 
    ctime_grid(1+(cc-1)*2,1:length(pcast)) = time_down; 

    % upcasts
    t1_grid(cc*2,1:length(pcast)) = t1up; 
    s1_grid(cc*2,1:length(pcast)) = s1up; 
    cp_grid(cc*2,1:length(pcast)) = pup; 
    ctime_grid(cc*2,1:length(pcast)) = time_up; 



    lons(cc,1:length(pcast)) = loncast; 
    lats(cc,1:length(pcast)) = latcast; 
    ctddates(cc,1:length(pcast)) = datecast; 
    %times(cc,1:length(pcast)) = timecast;     

    end

end

% what was the largest number of datapoints? 
disp(sprintf('longest array: %d',max(array_sizes)))

%% now plot 

figure()
subplot(2,2,1); hold on 
title('temp')
pcolor(ctddates,p_all,t1_all); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
xtickformat('dd HH:mm')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('thermal')); 
cb = colorbar;
ylabel(cb,'[deg C]','Rotation',270)
% clim([-0.2 0.2])

subplot(2,2,2); hold on 
title('sal')
pcolor(ctddates,p_all,s1_all); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('haline')); 
cb = colorbar;
ylabel(cb,'[psu]','Rotation',270)
clim([32 36])

subplot(2,2,3); hold on 
title('temp')
pcolor(ctime_grid,cp_grid,t1_grid); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
%xtickformat('dd HH:mm')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('thermal')); 
cb = colorbar;
ylabel(cb,'[deg C]','Rotation',270)
% clim([-0.2 0.2])

subplot(2,2,4); hold on 
title('sal')
pcolor(ctime_grid,cp_grid,s1_grid); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
%xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('haline')); 
cb = colorbar;
ylabel(cb,'[psu]','Rotation',270)
clim([32 36])




%% but for now just save designated values 

% station 1 

CTD1_1m = struct; 
CTD1_1m.lons = lons; CTD1_1m.lats = lats; CTD1_1m.times = times; CTD1_1m.dates = dates; 
CTD1_1m.t1 = t1_all; CTD1_1m.s1 = s1_all; CTD1_1m.p = p_all; 
CTD1_1m.t2 = t2_all; CTD1_1m.s2 = s2_all; CTD1_1m.depth = depth_all; 

save([out_dir 'CTD_1m_station' num2str(station)],'CTD1_1m'); 

%% save ctd depths to interpolate 

ctd_p = p_all; 


%% interpolate LADCP and 1m CTD onto 4m grid 

% check the resolution of each 

figure(); hold on

ladcp_dz = diff(p_grid(1,:)); 
ctd_dz = diff(cp_grid(1,:)); 


ladcphist = histogram(ladcp_dz,50);
ladcppdf = histcounts(ladcp_dz,50,'Normalization','pdf');

ctdhist = histogram(ctd_dz,50);
ctdpdf = histcounts(ctd_dz,50,'Normalization','pdf');

figure(); hold on 
ladcpbins = ladcphist.BinEdges + (ladcphist.BinWidth/2);
plot(ladcpbins(1:end-1), ladcppdf, 'b-')

ctdbins = ctdhist.BinEdges + (ctdhist.BinWidth/2);
plot(ctdbins(1:end-1), ctdpdf, 'r-')

%%
figure(); hold on
histogram(ladcp_dz,'Normalization','pdf','FaceColor','red','DisplayName',sprintf(' %d',mean(ladcp_dz,'all','omitnan')))
histogram(ctd_dz,'Normalization','pdf','FaceColor','blue','DisplayName',sprintf(' %d',mean(ctd_dz,'all','omitnan')))
legend()

%% 

% so lets interpolate onto the ladcp grid 

ladcp_z = p_grid(1,~isnan(p_grid(1,:)));
[nt,nz] = size(p_grid);
nz = length(ladcp_z);

p_gridded = repmat(ladcp_z,nt,1);
t1_gridded = nan(nt,nz);
s1_gridded = nan(nt,nz);

u_gridded = nan(nt,nz);
v_gridded = nan(nt,nz);

time_gridded = nan(nt,nz);

% interpolate each cast individually to the ladcp grid 
for cc = 1:nt
    % for 2d: vq = griddata()
    % p_gridded(cc,:) = interp1(ctd_p(cc,:),ctd_p(cc,:),ladcp_z);
    t1_gridded(cc,:) = interp1(cp_grid(cc,~isnan(cp_grid(cc,:))),t1_grid(cc,~isnan(cp_grid(cc,:))),ladcp_z);
    s1_gridded(cc,:) = interp1(cp_grid(cc,~isnan(cp_grid(cc,:))),s1_grid(cc,~isnan(cp_grid(cc,:))),ladcp_z);

    % also interpolate LADCP to this grid because of missing bottom data and likely inconsistencies 
    u_gridded(cc,:) = interp1(p_grid(cc,~isnan(p_grid(cc,:))),u_grid(cc,~isnan(p_grid(cc,:))),ladcp_z);
    v_gridded(cc,:) = interp1(p_grid(cc,~isnan(p_grid(cc,:))),v_grid(cc,~isnan(p_grid(cc,:))),ladcp_z);
    p_gridded(cc,:) = interp1(p_grid(cc,~isnan(p_grid(cc,:))),p_grid(cc,~isnan(p_grid(cc,:))),ladcp_z);
    time_on_grid(cc,:) = interp1(p_grid(cc,~isnan(p_grid(cc,:))),time_grid(cc,~isnan(p_grid(cc,:))),ladcp_z);

    % time_gridded(cc,:) = interp1(ctd_p(cc,~isnan(ctd_p(cc,:))),times(cc,~isnan(ctd_p(cc,:))),ladcp_z);
    
end

%time_dim = repmat(time,1,nz); 
%time_dim = datetime(time_gridded,'ConvertFrom','posixtime'); 
p_dim = repmat(ladcp_z,nt,1); 

figure(); hold on 
subplot(2,1,1); hold on
pcolor(time_on_grid,p_gridded,t1_gridded); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('thermal')); 
cb = colorbar;
ylabel(cb,'[degC]','Rotation',270)

subplot(2,1,2); hold on
pcolor(time_on_grid,p_gridded,s1_gridded); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('haline')); 
cb = colorbar;
ylabel(cb,'[psu]','Rotation',270)


figure(); hold on 
subplot(2,1,1); hold on
pcolor(time_on_grid,p_gridded,u_gridded); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('balance')); 
cb = colorbar;
ylabel(cb,'[m/s]','Rotation',270)
clim([-0.4 0.4])

subplot(2,1,2); hold on
pcolor(time_on_grid,p_gridded,v_gridded); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('balance')); 
cb = colorbar;
ylabel(cb,'[m/s]','Rotation',270)
clim([-0.4 0.4])


%% save that! cuz now you can do richardson calculations 



%% but also now make a least squares harmonic fit 
% and save that too 

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



%% lets do richardson calculation 

addpath /home/vboatwright/OneDrive/Documents/SIO/projects/gsw_matlab
addpath /home/vboatwright/OneDrive/Documents/SIO/projects/gsw_matlab/library
addpath /home/vboatwright/OneDrive/Documents/SIO/projects/gsw_matlab/thermodynamics_from_t

% first: change fields from practical salinity to absolute salinity 

SA_gridded = gsw_SA_from_SP(s1_gridded,p_gridded,lon_cast,lat_cast); 
CT_gridded = gsw_CT_from_t(SA_gridded,t1_gridded,p_gridded); 
p_ref = 0; 
rho_gridded = gsw_rho(SA_gridded,CT_gridded,p_ref); 

[nt,nz] = size(p_gridded);
tot_S2 = nan(nt,nz-1); 

for cc = 1:nt 
    % need to do individually for each column because they each start at a different point due to NaNs
    
    for ii = flip(2:nz)
        % going from the bottom to top  
        
        if isnan(u_grid(cc,ii-1)) == 0
            % do shear calculations only if the level above has a value 
            shear_u = (u_gridded(cc,ii-1)-u_gridded(cc,ii)) / (p_gridded(cc,ii-1) - p_gridded(cc,ii)); 
            shear_v = (v_gridded(cc,ii-1)-v_gridded(cc,ii)) / (p_gridded(cc,ii-1) - p_gridded(cc,ii));
            
            tot_S2(cc,ii-1) = shear_u^2 + shear_v^2;
        end 
    end
end


% Ri = N2 / S2 
% mixing at high ri-1 
% N2 = -g/rho d rho / dz 
% [N2, p_mid] = gsw_Nsquared(SA,CT,p,{lat})

[N2, p_mid] = gsw_Nsquared(SA_gridded',CT_gridded',p_gridded',lat_cast);

inv_Ri = tot_S2'./N2;  


%%
figure(); hold on 
subplot(2,1,1); hold on
pcolor(time_on_grid,p_gridded,rho_gridded); hold on
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('dense')); 
cb = colorbar;
ylabel(cb,'[kg/m3]','Rotation',270)
xticks(time_on_grid(1:2:end,1))
xticklabels(datestr(time_on_grid(1:2:end,1),'mm/dd HH:MM'))


timenums = datenum(time_on_grid);
clevels = [1];

subplot(2,1,2); hold on
pcolor(timenums(:,1:end-1),p_gridded(:,1:end-1),inv_Ri')
contour(timenums(:,1:end-1),p_gridded(:,1:end-1),inv_Ri',clevels,'Color','k')
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
colormap(gca,cmocean('matter')); 
cb = colorbar;
ylabel(cb,'Ri^{-1} #','Rotation',270)
set(gca,'ColorScale','log')
xticks(timenums(1:2:end,1))
xticklabels(datestr(time_on_grid(1:2:end,1),'mm/dd HH:MM'))

%%
% 
% for cc = 1:nt
%     for zz = 1:nz
% 
%         if isnan(u_gridded(cc,ii+1)) == 0
%             % do shear calculations only if the level below has a value 
%             shear_u = (u_gridded(cc,ii+1)-u_gridded(cc,ii)) / (p_gridded(cc,ii+1) - p_gridded(cc,ii)); 
%             shear_v = (v_gridded(cc,ii+1)-v_gridded(cc,ii)) / (p_gridded(cc,ii+1) - p_gridded(cc,ii));
% 
%             tot_S2 = shear_u^2 + shear_v^2; 
%         tot_S2 = 
%         [N2, p_mid] = gsw_Nsquared(SA_gridded,CT_gridded,p_gridded,lat_cast);
%         % manually: 
%         % n2 = N2[1:,:]-N2[0:-1,:])/(t[1:,None]-t[0:-1,None]) # dN2_dt[mid_t, mid_z]
%         ri_inv =  
%     end
% end


%%  station 2

CTD2 = struct; 
CTD2.lons = lons; CTD2.lats = lats; CTD2.times = times; CTD2.dates = dates; 
CTD2.t1 = t1_all; CTD2.s1 = s1_all; CTD2.p = p_all; 
CTD2.t2 = t2_all; CTD2.s2 = s2_all; CTD2.depth = depth_all; 

save([out_dir 'CTD_station' num2str(station)],'CTD2'); 
