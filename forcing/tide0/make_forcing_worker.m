function make_forcing_worker(grid_dir, tide_dir, out_dir, date_string)

% To test this function alone comment out the function declaration above
% and uncomment these four lines.
% grid_dir = '/Users/pm8/Documents/LO_data/grids/cas6';
% tide_dir = '/Users/pm8/Documents/LO_data/tide';
% out_dir = '/Users/pm8/Documents/LO_output/cas6_v3/f2019.07.04/tide0';
% date_string = '2019.07.04';

% make_forcing_worker.m
%
% ****************** for tide ***********************
%
% This makes the tidal NetCDF forcing file for a ROMS 3.# simulation, and
% NOTE that you need the TPXO model files (tmd_toolbox/)
%
% This version 2017.11.02 includes an increase to the amplitudes.
%
% 2021.03.23 LO Version:
% Pass directory names as strings with no trailing /.
% The main difference between this and the LiveOcean/forcing/tide2 original
% is that I converted all the NetCDF operations to use the matlab-native
% versions. This leaves me uneasy about default packing of things while
% reading and writing.


%% tide-specific code

% define needed files
gridfile = [grid_dir,'/grid.nc'];
t_dir = [tide_dir,'/TPXO/'];
ncfile_out = [out_dir,'/tides.nc'];  % tide forcing file name

if exist(ncfile_out, 'file')==2
    delete(ncfile_out);
end

%% parse the time in date_string for the nodal corrections
yrs = date_string(1:4);
mos = date_string(6:7);
dys = date_string(9:10);
yr = str2double(yrs);
mo = str2double(mos);
dy = str2double(dys);
tref_datenum = datenum(yr,mo,dy);

% get sizes and initialize the NetCDF file
% specify the constituents to use
c_data =  ['m2  ';'s2  ';'k1  ';'o1  '; ...
    'n2  ';'p1  ';'k2  ';'q1  '; ...
    '2n2 ';'mu2 ';'nu2 ';'l2  '; ...
    't2  ';'j1  ';'no1 ';'oo1 '; ...
    'rho1';'mf  ';'mm  ';'ssa ';'m4  '];
if 0
    np = size(c_data,1);
else
    % limit number of constituents, this is what is in TPXO7.1
    np = 8;
end

% note the transpose after reading
lon_rho = ncread(gridfile,'lon_rho')';
lat_rho = ncread(gridfile,'lat_rho')';
mask_rho = ncread(gridfile,'mask_rho')';

[ny,nx] = size(lon_rho);

% get constituent field and interpolate to grid

% get extracted amplitude and phases
t_file = [t_dir,'Extractions/tpxo7p2_180.mat'];
load(t_file);
x = Z.lon; y = Z.lat; mask = Z.mask;
% get set up to find constituent and nodal correction info
path(path,[t_dir,'/tmd_toolbox']);
% specify the tidal model to use
mod_name = 'tpxo7.2';
model = [t_dir,'/DATA7p2/Model_',mod_name];
for ii = 1:np
    cons = c_data(ii,:);
    cons_nb = deblank(cons);
    [junk,junk,ph,om,junk,junk] = ...
        constit(cons);
    disp([' Working on ',cons_nb,' period = ', ...
        num2str(2*pi/(3600*om)),' hours']);
    eval(['Eamp = E_',cons_nb,'.amp;']);
    eval(['Ephase = E_',cons_nb,'.phase;']);
    eval(['Cangle = C_',cons_nb,'.uincl;']);
    eval(['Cphase = C_',cons_nb,'.uphase;']);
    eval(['Cmax = C_',cons_nb,'.umaj;']);
    eval(['Cmin = C_',cons_nb,'.umin;']);
    % get nodal corrections centered on the middle
    % of the run time period, relative to day 48622 mjd (1/1/1992)
    trel = tref_datenum - datenum(1992,1,1);
    [pu,pf] = nodal(trel + 48622,cons);
    % interpolate to the model grid
    EC_list = {'Eamp','Ephase','Cangle','Cphase','Cmax','Cmin'};
    [X,Y] = meshgrid(x,y);
    for jj = 1:length(EC_list)
        eval(['this_in = ',EC_list{jj},';']);
        % first do a standard interpolation with NaN's
        this_out = NaN*lon_rho;
        this_in(~mask) = NaN;
        this_out = interp2(X,Y,this_in,lon_rho,lat_rho);
        % then go through row by row and extrapolate E-W
        nn = size(lon_rho,1);
        for kk = 1:nn
            this_row = this_out(kk,:);
            this_lon = lon_rho(kk,:);
            this_mask = ~isnan(this_row);
            if(length(this_row(this_mask))==1)
                % only 1 good element, but need 2 for interp1
                this_out_ex(kk,:) = this_row(this_mask);
            elseif(~isempty(this_row(this_mask)) && ...
                    length(this_row(this_mask))>1) % DAS added
                this_out_ex(kk,:) = interp1(this_lon(this_mask), ...
                    this_row(this_mask),this_lon,'nearest','extrap');
            else
                this_out_ex(kk,:) = this_out_ex(kk-1,:);
            end
        end
        eval([EC_list{jj},' = this_out_ex;']);
    end
    % and make final adjustments before writing to arrays
    disp([cons_nb, ': pf before = ',num2str(pf)])
    % PM Edit: diurnals
    if cons_nb == 'o1'
        pf = pf*1.21*1.087;
        phase_shift = -10.0; % deg
    elseif cons_nb == 'k1'
        pf = pf*1.21*1.11;
        phase_shift = -18.0; % deg
    elseif cons_nb == 'p1'
        pf = pf*1.21;
        phase_shift = -18.0; % deg
    elseif cons_nb == 'q1'
        pf = pf*1.21;
        phase_shift = -18.0; % deg
        % PM Edit: semidiurnals
    elseif cons_nb == 'm2'
        pf = pf*1.17*1.075;
        phase_shift = -25.0; % deg
    elseif cons_nb == 's2'
        pf = pf*1.261*1.13;
        phase_shift = -35.0; % deg
    elseif cons_nb == 'n2'
        pf = pf*1.196*1.11;
        phase_shift = -23.0; % deg
    elseif cons_nb == 'k2'
        pf = pf*1.2*1.11;
        phase_shift = -23.0; % deg
    end
    
    disp([cons_nb, ': pf after = ',num2str(pf)])
    tide_period(ii) = 2*pi/(3600*om); % hours
    tide_Eamp(ii,:,:) = pf*Eamp; % m
    tide_Ephase(ii,:,:) = Ephase - 180*ph/pi - 180*pu/pi;% + phase_shift; % deg
    tide_Cangle(ii,:,:) = Cangle; % deg
    tide_Cphase(ii,:,:) = Cphase - 180*ph/pi - 180*pu/pi;% + phase_shift; % deg
    tide_Cmax(ii,:,:) = pf*Cmax/100; % m s-1
    tide_Cmin(ii,:,:) = pf*Cmin/100; % m s-1
    
end

% make sure tide_period is a column vector
tide_period = tide_period(:);

disp(['*** creating ',ncfile_out,' ***']);

nccreate(ncfile_out, 'tide_period', 'Dimensions', {'tide_period', np});
ncwrite(ncfile_out, 'tide_period', tide_period);

nccreate(ncfile_out, 'tide_Eamp', 'Dimensions', {'xi_rho', nx, 'eta_rho', ny,'tide_period', np});
ncwrite(ncfile_out, 'tide_Eamp', permute(tide_Eamp,[3,2,1]));

nccreate(ncfile_out, 'tide_Ephase', 'Dimensions', {'xi_rho', nx, 'eta_rho', ny,'tide_period', np});
ncwrite(ncfile_out, 'tide_Ephase', permute(tide_Ephase,[3,2,1]));

nccreate(ncfile_out, 'tide_Cphase', 'Dimensions', {'xi_rho', nx, 'eta_rho', ny,'tide_period', np});
ncwrite(ncfile_out, 'tide_Cphase', permute(tide_Cphase,[3,2,1]));

nccreate(ncfile_out, 'tide_Cangle', 'Dimensions', {'xi_rho', nx, 'eta_rho', ny,'tide_period', np});
ncwrite(ncfile_out, 'tide_Cangle', permute(tide_Cangle,[3,2,1]));

nccreate(ncfile_out, 'tide_Cmax', 'Dimensions', {'xi_rho', nx, 'eta_rho', ny,'tide_period', np});
ncwrite(ncfile_out, 'tide_Cmax', permute(tide_Cmax,[3,2,1]));

nccreate(ncfile_out, 'tide_Cmin', 'Dimensions', {'xi_rho', nx, 'eta_rho', ny,'tide_period', np});
ncwrite(ncfile_out, 'tide_Cmin', permute(tide_Cmin,[3,2,1]));

% ncdisp(ncfile_out);



