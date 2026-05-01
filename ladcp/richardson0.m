%% calculating Ri 

% for Richardson #, need stratification and total shear 


% get LADCP & CTD data 
addpath '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/st_lucia_analysis/' 

figure_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/figures/ladcp_figures';

data_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/processed_ship/'; 
fn = 'LADCP_station1.mat'; 

file = [data_dir fn]; 
load(file) 

[nt,nz] = size(LADCP1.u);
u = LADCP1.u; v = LADCP1.v; ladcp_p = LADCP1.p; 
ladcp_time = LADCP1.time;

ladcptime_dim = repmat(ladcp_time,1,nz); 

fig1 = figure();
set(gcf, 'Color', 'w', 'Position', [100, 100, 900, 600]); % white background + larger size
clims = [-0.25, 0.25];

% u velocity 
subplot(2,1,1); hold on
pcolor(ladcptime_dim, ladcp_p, u); shading flat
colormap(cmocean('balance')); 
cb1 = colorbar;
ylabel(cb1, '[m/s]', 'Rotation', 270, 'FontSize', 11)
xlabel(''); 
ylabel('\bfDepth\rm [dbar]', 'FontSize', 12)
title({'\bfSLE Station 1', 'u velocity'}, 'FontSize', 13)
set(gca, 'YDir', 'reverse', 'FontSize', 11, 'LineWidth', 1)
xtickformat('MM-dd HH:mm')
caxis(clims)

% v velocity 
subplot(2,1,2); hold on
pcolor(ladcptime_dim, ladcp_p, v); shading flat
colormap(cmocean('balance')); 
cb2 = colorbar;
ylabel(cb2, '[m/s]', 'Rotation', 270, 'FontSize', 11)
xlabel('\bfDate\rm', 'FontSize', 12)
ylabel('\bfDepth\rm [dbar]', 'FontSize', 12)
title('v velocity', 'FontSize', 13)
set(gca, 'YDir', 'reverse', 'FontSize', 11, 'LineWidth', 1)
xtickformat('MM-dd HH:mm')
caxis(clims)

% can link axes for consistent zoom/pan
linkaxes(findall(gcf, 'Type', 'axes'), 'x')

%% based on wavenumber analysis, need to interpolate to 20m to calculate shear

% can we do a rolling smoother or do we need to straight up coarsen? 
% lets try both 


[nc,nz] = size(u);
rolling_shear_u = nan(nc,nz);
rolling_shear_v = nan(nc,nz);
nz_20m = 5; % roughly 5 datapoints in 20m 

bottom_depth = max(ladcp_p,[],'all'); 
shear_grid = 0:20:bottom_depth; 
interp_nz = length(shear_grid);

interp_u = nan(nc,interp_nz);
interp_v = nan(nc,interp_nz);

interp_shear_u = nan(nc,interp_nz);
interp_shear_v = nan(nc,interp_nz);



for cc = 1:nc 
    % need to do individually for each column because they each start at a different point due to NaNs
    
    % rolling method 
    deepest = sum(~isnan(ladcp_p(cc,:)));
    for ii = 3:deepest-3

        % going from the top to bottom - calculate rolling
        % 5 datapoints in 20m --> take 2 on each side 
        
        u_i = mean(u(cc,ii-2:ii+2),'omitnan'); 
        p_i = mean(ladcp_p(cc,ii-2:ii+2),'omitnan');
        u_ii = mean(u(cc,ii-1:ii+3),'omitnan');
        p_ii = mean(ladcp_p(cc,ii-1:ii+3),'omitnan');

        rolling_shear_u(cc,ii) = ( u_ii-u_i) / (p_ii-p_i ); 
        
    end

    % coarsening method 
    for ii = 1:interp_nz
        depth_on_grid = shear_grid(ii);
        [~,depth_index] = min(abs(ladcp_p(1,:)-depth_on_grid));

        if ii == 1
            % we are at the surface - only use top 2 in the mean 
            interp_u(cc,ii) = mean(u(cc,depth_index:depth_index+2),'omitnan');

        elseif ii == interp_nz
            % we are at the bottom - only use the bottom available in the mean
            interp_u(cc,ii) = mean(u(cc,depth_index-2:end),'omitnan');

        else 
            interp_u(cc,ii) = mean(u(cc,depth_index-2:depth_index+2),'omitnan');

        end

    end 

    % now, go through that loop but to calculate shear 
    for ii = 1:interp_nz-1
        interp_shear_u(cc,ii) = (interp_u(cc,ii+1)-interp_u(cc,ii)) / (shear_grid(ii+1) - shear_grid(ii)); 

    end
end


% old: 
% if isnan(u(cc,ii+1)) == 0
%     % do shear calculations only if the level below has a value 
%     %shear_u(cc,ii) = (u(cc,ii+1)-u(cc,ii)) / (p(cc,ii+1) - p(cc,ii)); 
%     %shear_v(cc,ii) = (u(cc,ii+1)-u(cc,ii)) / (p(cc,ii+1) - p(cc,ii));
% 
% 
% end 

%% now plot

[x,y] = meshgrid(ladcp_time,shear_grid);

fig1=figure()
subplot(3,1,1); hold on 
[t,s] = title('Station 1','interpolated u velocity');
pcolor(x',y',interp_u); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'u velocity [m/s]','Rotation',270)
clim([-0.25 0.25])

subplot(3,1,2); hold on 
subtitle('interpolated u shear')
pcolor(x',y',interp_shear_u); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'du/dz [1/s]','Rotation',270)
clim([-0.01 0.01])


subplot(3,1,3); hold on 
subtitle('rolling u shear')
pcolor(ladcptime_dim, ladcp_p,rolling_shear_u); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'du/dz [1/s]','Rotation',270)
clim([-0.01 0.01])


%%


subplot(2,1,2); hold on 
subtitle('v shear')
pcolor(time_dim,p,shear_v); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'dv/dz [1/s]','Rotation',270)
clim([-0.01 0.01])

saveas(fig1,sprintf('%s/postation1_ladcp_uv_shear.png',figure_dir),'png')


Z_coarse_avg = arrayfun(@(i,j) mean(mean(Z_fine(i:i+block_size-1, j:j+block_size-1))), ...
                                 1:block_size:size(Z_fine,1), 1:block_size:size(Z_fine,2));





%%
% and CTD data 
fn = 'CTD_station1.mat'; 

file = [data_dir fn]; 
load(file) 

[nt,nz] = size(CTD1.t1);
t1 = CTD1.t1; s1 = CTD1.s1; ctd_p = CTD1.p; 
t2 = CTD1.t2; s2 = CTD1.s2; depth = CTD1.depth; 
ctd_dates = CTD1.dates;

% temp
subplot(2,2,1); hold on
pcolor(ctd_dates, ctd_p, t1); shading flat
colormap(gca,cmocean('thermal')); 
cb1 = colorbar;
ylabel(cb1, 'degC', 'Rotation', 270, 'FontSize', 11)
xlabel('temperature 1'); 
ylabel('\bfDepth\rm [dbar]', 'FontSize', 12)
title({'\bfSLE Station 1', 'temperature'}, 'FontSize', 13)
set(gca, 'YDir', 'reverse', 'FontSize', 11, 'LineWidth', 1)
xtickformat('MM-dd HH:mm')

% salt 
subplot(2,2,2); hold on
pcolor(ctd_dates, ctd_p, s1); shading flat
colormap(gca,cmocean('haline')); 
cb2 = colorbar;
ylabel(cb2, 'psu', 'Rotation', 270, 'FontSize', 11)
xlabel('\bfDate\rm', 'FontSize', 12)
ylabel('\bfDepth\rm [dbar]', 'FontSize', 12)
title('salinity 1', 'FontSize', 13)
set(gca, 'YDir', 'reverse', 'FontSize', 11, 'LineWidth', 1)
xtickformat('MM-dd HH:mm')

% and also second ctd 

% temp
subplot(2,2,3); hold on
pcolor(ctd_dates, ctd_p, t2); shading flat
colormap(gca,cmocean('thermal')); 
cb1 = colorbar;
ylabel(cb1, 'degC', 'Rotation', 270, 'FontSize', 11)
xlabel('temperatyre 2'); 
ylabel('\bfDepth\rm [dbar]', 'FontSize', 12)
title({'\bfSLE Station 1', 'temperature'}, 'FontSize', 13)
set(gca, 'YDir', 'reverse', 'FontSize', 11, 'LineWidth', 1)
xtickformat('MM-dd HH:mm')

% salt 
subplot(2,2,4); hold on
pcolor(ctd_dates, ctd_p, s2); shading flat
colormap(gca,cmocean('haline')); 
cb2 = colorbar;
ylabel(cb2, 'psu', 'Rotation', 270, 'FontSize', 11)
xlabel('\bfDate\rm', 'FontSize', 12)
ylabel('\bfDepth\rm [dbar]', 'FontSize', 12)
title('salinity 2', 'FontSize', 13)
set(gca, 'YDir', 'reverse', 'FontSize', 11, 'LineWidth', 1)
xtickformat('MM-dd HH:mm')


% can link axes for consistent zoom/pan
linkaxes(findall(gcf, 'Type', 'axes'), 'x')

%% also wtf is happening with t2 


% temp
masked_t2 = t2;  
masked_t2(masked_t2 <= 0) = NaN; 

figure(); hold on
subplot(2,1,1)
pcolor(ctd_dates, p, masked_t2); shading flat
colormap(gca,cmocean('thermal')); 
cb1 = colorbar;
ylabel(cb1, 'degC', 'Rotation', 270, 'FontSize', 11)
xlabel('temperature '); 
ylabel('\bfDepth\rm [dbar]', 'FontSize', 12)
title({'\bfSLE Station 1', 'temperature above 0 deg C'}, 'FontSize', 13)
set(gca, 'YDir', 'reverse', 'FontSize', 11, 'LineWidth', 1)
xtickformat('MM-dd HH:mm')


% temp
masked_t2 = t2;  
masked_t2(masked_t2 >= 0) = NaN;  

subplot(2,1,2); hold on 
pcolor(ctd_dates, p, masked_t2); shading flat
colormap(gca,cmocean('thermal')); 
cb1 = colorbar;
ylabel(cb1, 'degC', 'Rotation', 270, 'FontSize', 11)
xlabel('temperature'); 
ylabel('\bfDepth\rm [dbar]', 'FontSize', 12)
title({'\bfSLE Station 1', 'temperature below 0 deg C'}, 'FontSize', 13)
set(gca, 'YDir', 'reverse', 'FontSize', 11, 'LineWidth', 1)
xtickformat('MM-dd HH:mm')

figure(); hold on 
subplot(2,1,1); hold on 
plot(p,t2)
ylabel('temperature'); 
xlabel('depth')
title('Temperature from casts across depths')
subplot(2,1,2); hold on 
plot(CTD1.times,t2)
%xticks(sort(CTD1.times))
%xticklabels(ladcp_time)
ylabel('temperature'); 
xlabel('time')
title('Temperature from casts over time')

%% interpolate ladcp and ctd to the same grid 

% use just t1 and s1 for now bc messed up t2/s2 





%% total shear for Richardson # 


[nc,nz] = size(u);
shear_u = nan(nc,nz);
shear_v = nan(nc,nz);
shear_squared = nan(nc,nz);

for cc = 1:nc 
    % need to do individually for each column because they each start at a different point due to NaNs
    
    for ii = 1:nz-1
        % going from the top to bottom  
        
        if isnan(u(cc,ii+1)) == 0
            % do shear calculations only if the level below has a value 
            shear_u(cc,ii) = (u(cc,ii+1)-u(cc,ii)) / (p(cc,ii+1) - p(cc,ii)); 
            shear_v(cc,ii) = (u(cc,ii+1)-u(cc,ii)) / (p(cc,ii+1) - p(cc,ii));
            shear_squared(cc,ii) = (shear_u(cc,ii)^2) + (shear_v(cc,ii) )^2 ; 

        end 
    end
end


% now plot

fig1=figure()
subplot(3,1,1); hold on 
[t,s] = title('Station 1','u shear');
pcolor(time_dim,p,shear_u); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'du/dz [1/s]','Rotation',270)
clim([-0.01 0.01])

subplot(3,1,2); hold on 
subtitle('v shear')
pcolor(time_dim,p,shear_v); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'dv/dz [1/s]','Rotation',270)
clim([-0.01 0.01])


subplot(3,1,3); hold on 
subtitle('total shear')
pcolor(time_dim,p,shear_tot); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'total shear [1/s]','Rotation',270)
%clim([0 0.02])



%% next calculate N2 

% need TS! load CTD mat file 




