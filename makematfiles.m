% script for concatenating and processing u/v velocities from SLE stations 

% data_dir = /change this to kerstin's harddrive ; 
data_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/processed_ship/'; 
out_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/processed_ship/';

station = 2; 


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


%% process and concatenate 

% for station 1: 18 casts 
if station == 1 
    casts = 1:18; 
elseif station == 2
    casts = 1:36; 
end
    
numcasts = length(casts);

% preallocate based on cast 1 

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

lon = nan(numcasts,1); 
lat = nan(numcasts,1); 
time = NaT(numcasts,1);


for cc = casts
 
    if cc == 21 || cc ==26 
        % skip these 
        % leave as nan 

    else

    file = sprintf('%02d.mat',cc);
    fn = [ladcp_dir file]; 
    S = load(fn); 
    
    % load variables
    u_cast = S.dr.u; 
    v_cast = S.dr.v; 

    ushear = S.dr.u_shear_method; 
    vshear = S.dr.v_shear_method; 
    wshear = S.dr.w_shear_method; 
    
    
    lon_cast = S.dr.lon; 
    lat_cast = S.dr.lat; 
    
    p_cast = S.dr.p; 
    time_cast = datetime(S.dr.date); 
     
    %disp(length(p_cast))
    
    % check the uv ctd 
    u_ctd_cast = S.dr.uctd; 
    v_ctd_cast = S.dr.vctd; 
    w_ctd_cast = S.dr.wctd; 
    
    % ship info 

    slon = S.dr.shiplon; 
    slat = S.dr.shiplat; 

    disp(length(u_ctd_cast))

    % save T, S - though bad data 
    t_ctd = S.dr.ctd_t;
    s_ctd = S.dr.ctd_s;

    % concatenate 

    lon(cc) = lon_cast; 
    lat(cc) = lat_cast;
    time(cc) = time_cast;

    u_all(cc,1:length(p_cast)) = u_cast; 
    v_all(cc,1:length(p_cast)) = v_cast; 
    p_all(cc,1:length(p_cast)) = p_cast; 

    u_shearmethod(cc,1:length(p_cast)) = ushear; 
    v_shearmethod(cc,1:length(p_cast)) = vshear; 
    w_shearmethod(cc,1:length(p_cast)) = wshear; 


    t_all(cc,1:length(p_cast)) = t_ctd; 
    s_all(cc,1:length(p_cast)) = s_ctd; 
    
    u_ctd(cc,1:length(u_ctd_cast)) = u_ctd_cast; 
    v_ctd(cc,1:length(u_ctd_cast)) = u_ctd_cast;
    w_ctd(cc,1:length(u_ctd_cast)) = u_ctd_cast;

    ship_lon(cc,1:length(u_ctd_cast)) = slon;    
    ship_lat(cc,1:length(u_ctd_cast)) = slat;

    end

end

%% now plot 

cast_dim = repmat(casts,zlevels,1)'; 
time_dim = repmat(time,1,zlevels); 

figure()
subplot(2,1,1); hold on 
title('u velocity')
pcolor(time_dim,p_all,u_all); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
xtickformat('dd HH:mm')
set(gca, 'YDir', 'reverse')
cb = colorbar;
ylabel(cb,'u [m/s]','Rotation',270)
clim([-0.2 0.2])

subplot(2,1,2); hold on 
title('v velocity')
pcolor(time_dim,p_all,v_all); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
cb = colorbar;
ylabel(cb,'v [m/s]','Rotation',270)
clim([-0.2 0.2])


%% calculate shear 

[nc,nz] = size(u_all);
shear_u = nan(nc,nz);
shear_v = nan(nc,nz);

disp(size(shear_u))
for cc = casts 
    % need to do individually for each column because they each start at a different point due to NaNs
    
    for ii = 1:nz-1
        % going from the top to bottom  
        
        if isnan(u_all(cc,ii+1)) == 0
            % do shear calculations only if the level below has a value 
            shear_u(cc,ii) = (u_all(cc,ii+1)-u_all(cc,ii)) / (p_all(cc,ii+1) - p_all(cc,ii)); 
            shear_v(cc,ii) = (v_all(cc,ii+1)-v_all(cc,ii)) / (p_all(cc,ii+1) - p_all(cc,ii));

        end 
    end
end


%% now plot shear

figure()
subplot(2,1,1); hold on 
title('u shear')
pcolor(cast_dim,p_all,shear_u); 
shading flat 
xlabel('Cast #'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cb = colorbar;
ylabel(cb,'du/dz [1/s]','Rotation',270)
clim([-0.0125 0.0125])

subplot(2,1,2); hold on 
title('v shear')
pcolor(cast_dim,p_all,shear_v); 
shading flat 
xlabel('Cast #'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cb = colorbar;
ylabel(cb,'dv/dz [1/s]','Rotation',270)
clim([-0.0125 0.0125])


%% but for now just save designated values 

% station 1 

LADCP1 = struct; 
LADCP1.lon = lon; LADCP1.lat = lat; LADCP1.time = time; 
LADCP1.u = u_all; LADCP1.v = v_all; LADCP1.p = p_all; 
LADCP1.u_shearmethod = u_shearmethod; LADCP1.v_shearmethod = v_shearmethod; LADCP1.w_shearmethod = w_shearmethod; 
LADCP1.tctd = t_ctd; LADCP1.sctd = s_ctd; 
LADCP1.u_ctd = u_ctd; LADCP1.v_ctd = v_ctd; LADCP1.w_ctd = w_ctd; 
LADCP1.ship_lon = ship_lon; LADCP1.ship_lat = ship_lat; 

save([out_dir 'LADCP_station' num2str(station)],'LADCP1'); 

%% station 2 


LADCP2 = struct; 
LADCP2.lon = lon; LADCP2.lat = lat; LADCP2.time = time; 
LADCP2.u = u_all; LADCP2.v = v_all; LADCP2.p = p_all; 
LADCP2.u_shearmethod = u_shearmethod; LADCP2.v_shearmethod = v_shearmethod; LADCP2.w_shearmethod = w_shearmethod; 
LADCP2.tctd = t_ctd; LADCP2.sctd = s_ctd; 
LADCP2.u_ctd = u_ctd; LADCP2.v_ctd = v_ctd; LADCP2.w_ctd = w_ctd; 
LADCP2.ship_lon = ship_lon; LADCP2.ship_lat = ship_lat; 

save([out_dir 'LADCP_station' num2str(station)],'LADCP2'); 


%% now concatenate temp and salinity from .mat files 

ctd_dir = [data_dir 'ctd_mat/station' num2str(station) '/']; 

station = 2;
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

t1 = S.datad.t1; 
s1 = S.datad.s1; 
t2 = S.datad.t2; 
s2 = S.datad.s2; 

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

figure()
plot(t1,depth); hold on 
plot(s1,depth)

%% process and concatenate 

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

t1_all = nan(numcasts,zlevels); 
s1_all = nan(numcasts,zlevels); 
p_all = nan(numcasts,zlevels); 

t2_all = nan(numcasts,zlevels); 
s2_all = nan(numcasts,zlevels); 
depth_all = nan(numcasts,zlevels); 

lons = nan(numcasts,zlevels); 
lats = nan(numcasts,zlevels); 
dates = NaT(numcasts,zlevels); 
times = nan(numcasts,zlevels); 

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
    
    t1cast = S.datad.t1; 
    s1cast = S.datad.s1; 
    t2cast = S.datad.t2; 
    s2cast = S.datad.s2; 
    
    theta1cast = S.datad.theta1; 
    sigma1cast = S.datad.sigma1; 
    theta2cast = S.datad.theta2; 
    sigma2cast = S.datad.sigma2; 
    
    loncast = S.datad.lon; 
    latcast = S.datad.lat; 
    
    depthcast = S.datad.depth; 
    pcast = S.datad.p; 
    datecast = datetime(S.datad.datenum,'ConvertFrom','datenum'); 
    timecast = S.datad.time; 

    array_sizes = cat(1,array_sizes,length(pcast)); 

    % concatenate 

    t1_all(cc,1:length(pcast)) = t1cast; 
    s1_all(cc,1:length(pcast)) = s1cast; 
    p_all(cc,1:length(pcast)) = pcast; 
    
    t2_all(cc,1:length(pcast)) = t2cast; 
    s2_all(cc,1:length(pcast)) = s2cast; 
    depth_all(cc,1:length(pcast)) = depthcast; 

    lons(cc,1:length(pcast)) = loncast; 
    lats(cc,1:length(pcast)) = latcast; 
    dates(cc,1:length(pcast)) = datecast; 
    times(cc,1:length(pcast)) = timecast;     

    end

end

%% now plot 

figure()
subplot(2,1,1); hold on 
title('temp')
pcolor(dates,p_all,t1_all); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
xtickformat('dd HH:mm')
set(gca, 'YDir', 'reverse')
cb = colorbar;
ylabel(cb,'[deg C]','Rotation',270)
% clim([-0.2 0.2])

subplot(2,1,2); hold on 
title('sal')
pcolor(dates,p_all,s1_all); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
cb = colorbar;
ylabel(cb,'[psu]','Rotation',270)
clim([32 36])


%% but for now just save designated values 

% station 1 

CTD1 = struct; 
CTD1.lons = lons; CTD1.lats = lats; CTD1.times = times; CTD1.dates = dates; 
CTD1.t1 = t1_all; CTD1.s1 = s1_all; CTD1.p = p_all; 
CTD1.t2 = t2_all; CTD1.s2 = s2_all; CTD1.depth = depth_all; 

save([out_dir 'CTD_station' num2str(station)],'CTD1'); 

%%  station 2

CTD2 = struct; 
CTD2.lons = lons; CTD2.lats = lats; CTD2.times = times; CTD2.dates = dates; 
CTD2.t1 = t1_all; CTD2.s1 = s1_all; CTD2.p = p_all; 
CTD2.t2 = t2_all; CTD2.s2 = s2_all; CTD2.depth = depth_all; 

save([out_dir 'CTD_station' num2str(station)],'CTD2'); 
