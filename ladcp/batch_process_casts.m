%batch_process_casts.m
%
% process one or more LADCP casts from the SW Indian Ocean
% using version IX_4 of thurnherr/visbeck software
%
% jen, dec 07


% set some cruise specific parameters

cruise = 'HIPPO';
stn='3000';
jrootdir=['/users/jen/projects/swirm/Ladcp/data/'];
yearbase=2007; 

% which ctd profiles go with which ladcp files?
if strcmp(stn,'2001')
    ctdlist='dodo2';
elseif strcmp(stn,'3000')
    ctdlist=['dodo3';'dodo4';'dodo5';'dodo6'];
elseif strcmp(stn,'4001')
    ctdlist=['dodo7';'dodo8';'dodo9'];
end

for icast=1:size(ctdlist,1)
disp(['PROCESSING STATION ' stn ' and CTD Cast # ' int2str(icast)])
    
% starty by setting a few things to default state
clear f d dr p ps;			% blank slate
f.ladcpdo = ' ';			% required by [m/default.m]
default;				% load default parameters
p = setdefv(p,'checkpoints',[]);	% disable checkpointing by default
close all;
more off;



%time range of a particular cast (many in same file)
% get times from ctd file
ctdname=ctdlist(icast,:);
load([jrootdir 'ctd/' ctdname '.mat'])
min_ctd_time=julian(yearbase,1,0,0) + min(CTD.yday(:));
max_ctd_time=julian(yearbase,1,0,0) + max(CTD.yday(:));
p.time_start=gregoria(min_ctd_time-4/1440);
p.time_end=gregoria(max_ctd_time+4/1440);
disp(['CTD start time = ']); disp(p.time_start)
disp(['CTD stop time = ']); disp(p.time_end)
clear CTD

% where to find and save  files
f.checkpoints = sprintf([jrootdir 'checkpoints/'  cruise  stn]);
f.res         = sprintf([jrootdir 'processed/' cruise stn '_' int2str(icast)]);
f.ladcpdo     = sprintf([jrootdir 'raw/' cruise  stn '.000']);

% load ctd time series
f.ctd 		  = sprintf([jrootdir 'ctd/' ctdname '_clean_lo.cnv']);
f.ctd_header_lines=0;
f.ctd_fields_per_line=4;
f.ctd_time_field=1;
f.ctd_pressure_field=2;
f.ctd_temperature_field=3;
f.ctd_salinity_field=4;
f.ctd_time_base=0;

% add shipboard adcp data, 150 for now
f.sadcp=['/users/jen/projects/swirm/adcp_150/sadcp_150.mat'];
p.sadcp_dtok=0;

% add ship gps data from pcode
f.nav=['users/jen/projects/swirm/nav/2007dec01_02.txt'];
f.nav_header_lines=0;
f.nav_fields_per_line=3;
f.nav_time_field=1;
f.nav_lat_field=2;
f.nav_lon_field=3;
f.nav_time_base=1;


% add lat/lon if we're not getting it from gps
%p.poss=[34 9.142 56 44.52];
%p.pose=[34 9.142 56 44.52];


% other bookeeping
p.cruise_id = cruise;
p.ladcp_station = stn; 
p.whoami = 'jen';
p.name = sprintf([p.cruise_id(1:3) ' ' p.cruise_id(5:end) ': #' p.ladcp_station ' dual']);
p.savecdf = 0;
p.saveplot_png = [1:4 9 11 13:14];
p.saveplot = [];
p.btrk_ts =  10; 
p = setdefv(p,'savemat',1);
%p.edit_mask_dn_bins = [1];
%p.down_sn = 8062;
%p.up_sn = 8028;
ps.smallfac = [1 0]; 


%clear pk; 

process_cast(stn,1,0,'',f,p,pk,ps)

end % icast loop

