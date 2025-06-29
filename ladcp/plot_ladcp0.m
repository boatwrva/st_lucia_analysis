%% plotting LADCP 

analysis_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/st_lucia_analysis';
addpath analysis_dir 

figure_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/figures/ladcp_figures';

data_dir = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/processed_ship/'; 
fn = 'LADCP_station1.mat'; 

file = [data_dir fn]; 
load(file) 

[nt,nz] = size(LADCP1.u);
u = LADCP1.u; v = LADCP1.v; p = LADCP1.p; 
time = LADCP1.time;

time_dim = repmat(time,1,nz); 

fig1 = figure();
set(gcf, 'Color', 'w', 'Position', [100, 100, 900, 600]); % white background + larger size
clims = [-0.25, 0.25];

% u velocity 
subplot(2,1,1); hold on
pcolor(time_dim, p, u); shading flat
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
pcolor(time_dim, p, v); shading flat
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


saveas(fig1,sprintf('%s/postation1_ladcp_uv.png',figure_dir),'png')



fig1 = figure(); set(gcf,'Color','white')
subplot(2,1,1); hold on 
[t,s] = title('SLE Station 1','u velocity');
pcolor(time_dim,p,u); shading flat
cmocean('balance'); cb = colorbar;
ylabel(cb,'u [m/s]','Rotation',270)
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
xtickformat('MM-dd HH:mm')
clim([-0.25 0.25])

subplot(2,1,2); hold on 
subtitle('v velocity')
pcolor(time_dim,p,v); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'v [m/s]','Rotation',270)
clim([-0.25 0.25])


%% shear 

[nc,nz] = size(u);
shear_u = nan(nc,nz);
shear_v = nan(nc,nz);

for cc = 1:nc 
    % need to do individually for each column because they each start at a different point due to NaNs
    
    for ii = 1:nz-1
        % going from the top to bottom  
        
        if isnan(u(cc,ii+1)) == 0
            % do shear calculations only if the level below has a value 
            shear_u(cc,ii) = (u(cc,ii+1)-u(cc,ii)) / (p(cc,ii+1) - p(cc,ii)); 
            shear_v(cc,ii) = (u(cc,ii+1)-u(cc,ii)) / (p(cc,ii+1) - p(cc,ii));

        end 
    end
end


% now plot

fig1=figure()
subplot(2,1,1); hold on 
[t,s] = title('Station 1','u shear');
pcolor(time_dim,p,shear_u); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'du/dz [1/s]','Rotation',270)
clim([-0.01 0.01])

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




%% station 2 

fn = 'LADCP_station2.mat'; 

file = [data_dir fn]; 
load(file) 

[nt,nz] = size(LADCP2.u);
u = LADCP2.u; v = LADCP2.v; p = LADCP2.p; 
time = LADCP2.time;

time_dim = repmat(time,1,nz); 

fig1 = figure();
set(gcf, 'Color', 'w', 'Position', [100, 100, 900, 600]); % white background + larger size
clims = [-0.25, 0.25];

% u velocity 
subplot(2,1,1); hold on
pcolor(time_dim, p, u); shading flat
colormap(cmocean('balance')); 
cb1 = colorbar;
ylabel(cb1, '[m/s]', 'Rotation', 270, 'FontSize', 11)
xlabel(''); 
ylabel('\bfDepth\rm [dbar]', 'FontSize', 12)
title({'\bfSLE Station 2', 'u velocity'}, 'FontSize', 13)
set(gca, 'YDir', 'reverse', 'FontSize', 11, 'LineWidth', 1)
xtickformat('MM-dd HH:mm')
caxis(clims)

% v velocity 
subplot(2,1,2); hold on
pcolor(time_dim, p, v); shading flat
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


saveas(fig1,sprintf('%s/postation2_ladcp_uv.png',figure_dir),'png')


fig2 = figure(); set(gcf,'Color','white')
subplot(2,1,1); hold on 
[t,s] = title('Station 2','u velocity');
pcolor(time_dim,p,u); shading flat
cmocean('balance'); cb = colorbar;
ylabel(cb,'u [m/s]','Rotation',270)
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
xtickformat('MM-dd HH:mm')
clim([-0.25 0.25])

subplot(2,1,2); hold on 
subtitle('v velocity')
pcolor(time_dim,p,v); 
shading flat 
xlabel('Date'); ylabel('Depth [dbar]')
xtickformat('MM-dd HH:mm')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'v [m/s]','Rotation',270)
clim([-0.25 0.25])


%% shear again 

[nc,nz] = size(u);
shear_u = nan(nc,nz);
shear_v = nan(nc,nz);

for cc = 1:nc 
    % need to do individually for each column because they each start at a different point due to NaNs
    
    for ii = 1:nz-1
        % going from the top to bottom  
        
        if isnan(u(cc,ii+1)) == 0
            % do shear calculations only if the level below has a value 
            shear_u(cc,ii) = (u(cc,ii+1)-u(cc,ii)) / (p(cc,ii+1) - p(cc,ii)); 
            shear_v(cc,ii) = (u(cc,ii+1)-u(cc,ii)) / (p(cc,ii+1) - p(cc,ii));

        end 
    end
end


% now plot

fig2=figure();
subplot(2,1,1); hold on 
[t,s] = title('Station 2','u shear');
pcolor(time_dim,p,shear_u); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'du/dz [1/s]','Rotation',270)
clim([-0.01 0.01])

subplot(2,1,2); hold on 
subtitle('v shear')
pcolor(time_dim,p,shear_v); 
shading flat 
xlabel('Time'); ylabel('Depth [dbar]')
set(gca, 'YDir', 'reverse')
cmocean('balance'); cb = colorbar;
ylabel(cb,'dv/dz [1/s]','Rotation',270)
clim([-0.01 0.01])


saveas(fig2,sprintf('%s/postation2_ladcp_uv_shear.png',figure_dir),'png')



%% least squares for harmonic 

