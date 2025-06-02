set(0, 'DefaultFigureColor', 'w');
path = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/bathymetry/';

filename = 'gebco_2024_n49.0_s30.0_w-130.0_e-110.0.nc';
%B = ncread(filename);
fn = '/home/vboatwright/OneDrive/Documents/SIO/projects/santalucia/data/bathymetry/gebco_2024_n49.0_s30.0_w-130.0_e-115.0.nc'; 

B.lat = ncread(fn,'lat');
B.lon = ncread(fn,'lon');
B.z = ncread(fn,'elevation');


%%
fig1 = figure(1);
clf
set_figure_LUCIA(fig1)
fig1.Units = 'normalized';

ftsize = 16;

long_west = -122.5; long_east = -120;
lat_south = 34; lat_north = 35.5; 

latlim0 = [lat_south lat_north];
lonlim0 = [long_west long_east];

[LONG, LAT] = meshgrid(B.lon, B.lat);
ELEV = B.z';

c_lims =[-4000 1];

m_proj('lambert','lat',[lat_south, lat_north],'lon',[long_west long_east], 'clo' , mean(lonlim0), 'par' , latlim0, 'rec' , 'off');
m_contourf(LONG, LAT, B.z',[-4500:500:-4001 -4000:250:0 ],'ShowText','on','Labelspacing',250)
 
cmapland=flipud(bone(30));
colormap(gca,[m_colmap('blues',60);cmapland(11:end,:)]);  
ccc=m_colmap('blues',60);
colormap(ccc(5:45,:))

caxis(gca,[-4000 0]);  
c1 = colorbar;
c1.Label.String = '[m]';

hold on
for ii=1:9
        m_plot(FinalStations{ii,3}, FinalStations{ii,2},'ro','Linewidth',10)
end


m_coast('patch',[.4 .4 .4],'edgecolor','none'); % this is the fast version
%m_gshhs_h('patch',[.7 .7 .7],'edgecolor','k'); % this is the slow detailed version of the land

m_grid('box','fancy','tickdir','in','yaxislocation','left',...
            'xaxislocation','bottom','xlabeldir','end','ticklen',.02);

m_ruler([.1 .25], .1 , 3,'fontsize',14,'Color','w','fontcolor','w')

