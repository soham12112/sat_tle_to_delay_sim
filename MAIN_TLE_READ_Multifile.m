%% TWO-LINE-ELEMENT READER / PLOTTER
%
% -------------------------------------------------------------------------
% DESCRIPTION
%
% Given multiple NORAD Two-Line-Element (TLE) files, this matlab code plots
% the orbits of the satellites as well as the Earth. This is meant to be a
% simple orbit plotting / visulization tool. 
% 
%% WORKSPACE SETUP

clc
clear all
close all
format longG

% Add path to the Earth plotting function. 
addpath([pwd, '/PlotEarth/PlotEarth']);

% Add path to the TLE data (these are just some examples with a focus on
% GNSS). 
addpath([pwd, '/TLE_Files']);

%% LOAD PHYSICAL CONSTANTS INTO THE GLOBAL WORKSPACE

physical_constants

global mu

%% SELECT THE TLE FILES TO PLOT

% All of GNSS. 
%filenames = {'gps-ops-clean','galileo-full','sbas-all-clean',...
 %   'beidou-clean-complete','glo-ops-clean'};

filenames = {'iridium_mine'};
% GPS + SBAS.
% filenames = {'gps-ops-clean','sbas-clean','qzss-full'};

% GLONASS. 
% filenames = {'glo-ops-clean'};

% Completed Galileo constellation. 
% filenames = {'galileo-full'};

% BeiDou
% filenames = {'beidou-05-25-2017'};

% GPS + Iridium
% filenames = {'gps-ops-clean','iridium'};

% GPS + SBAS GEOS. 
% filenames = {'gps-ops-clean','sbas-reduced-clean'};

% Transit. 
% filenames = {'nnss'};

% QZSS. 
% filenames = {'qzss-full'};

% Planet labs.
% filenames = {'rapideye', 'skybox', 'planet_mc'};

% Create colors for plotting. 
colors = lines(length(filenames));
% colors = hsv(length(filenames));

%% DOPPLER CALCULATION SETUP
% --- User Equipment (UE) Location ---
lat_ue_deg = 37.36;  % [deg] North
lon_ue_deg = -122.03; % [deg] West
alt_ue_m = 0;       % [m] Altitude above ellipsoid

% --- WGS84 Constants ---
R_e_wgs = 6378137.0;        % WGS84 Earth radius [m]
f_wgs = 1/298.257223563;  % WGS84 flattening
e_sq_wgs = 2*f_wgs - f_wgs^2; % WGS84 eccentricity squared
omega_E = 7.2921150e-5;   % Earth rotation rate [rad/s]
c_light = 299792458;      % Speed of light [m/s]

% --- Satellite Frequency (Assuming Iridium L-Band) ---
f_carrier = 1620e6; % [Hz]

% --- Convert UE Lat/Lon/Alt to ECEF ---
lat_ue_rad = deg2rad(lat_ue_deg);
lon_ue_rad = deg2rad(lon_ue_deg);

N = R_e_wgs / sqrt(1 - e_sq_wgs * sin(lat_ue_rad)^2);
R_ue_ecef = [ (N + alt_ue_m) * cos(lat_ue_rad) * cos(lon_ue_rad);
              (N + alt_ue_m) * cos(lat_ue_rad) * sin(lon_ue_rad);
              (N * (1 - e_sq_wgs) + alt_ue_m) * sin(lat_ue_rad) ];
          
% Velocity of the UE (fixed on Earth's surface) in ECEF
V_ue_ecef = [ -omega_E * R_ue_ecef(2);
               omega_E * R_ue_ecef(1);
               0 ];

%% DECODE TLE DATA

% Plot the Earth. 
% If you want a color Earth, use 'neompa', 'BlueMarble'.
% If you want a black and white Earth, use 'neomap', 'BlueMarble_bw'.
% A smaller sample step gives a finer resolution Earth.
h = plotearth('neomap', 'BlueMarble_bw', 'SampleStep', 2);

% Choose a date, by default this chooses the current date/time.
simStart = datenum(clock);
% simStart = datenum('Jan 09 2017 00:00:00');

% Compute sidereal time. 
GMST = utc2gmst(datevec(simStart)); % [rad]

% Import the TLE data. 
for k = 1:length(filenames)
    % Get the orbital elements. 
    [coe] = two_line_elem_conv(horzcat(filenames{k}, '.txt'), 'all');
    
    % Find latest epoch (all others can be run up from here).
    coeDateNums = datenum(coe.date);
    [val, ind] = min(coeDateNums);
    
    % Define the max time from simStart.
    a = coe.a(1);
    n = sqrt(mu/a^3);
    tFinal = 2*pi/n*1.1/3600/24; % This gives us just over 1 orbit. 
    
    % Create a time vector. 
    tSim = linspace(simStart, simStart + tFinal, 200);
    
    % Allocate space.
    RSave = NaN(length(tSim), 3, length(coeDateNums));
    
    % Run through all of the satellites in the TLE file 
    % and compute trajectories for plotting.
    for i = 1:length(coeDateNums) % For each satellite
        for j = 1:length(tSim) % For each time step
            % Get the orbit data. 
            a = coe.a(i);
            n = sqrt(mu / a^3);
            e = coe.e(i);
            inc = coe.i(i) * pi / 180;
            RAAN = coe.RAAN(i) * pi / 180;
            omega = coe.omega(i) * pi / 180;
            M = coe.M(i) * pi / 180 + ...
                n * (tSim(j) - coeDateNums(i)) * 24 * 3600; 
            
            % Adjust RAAN such that we are consisten with Earth's current
            % orientation. This is a conversion to Longitude of the
            % Ascending Node (LAN). 
            RAAN = RAAN - GMST;
            
            % Convert to ECI and save the data.
            [X,~] = COE2RV(a, e, inc, RAAN, omega, M);
            RSave(j,:,i) = X';
        end
    end
    
    % Plot the orbit.
    for i = 1:length(coeDateNums)
        colorI = k;
        plot3(RSave(:,1,i) / R_e, RSave(:,2,i) / R_e, RSave(:,3,i) / R_e,...
            'color', colors(colorI,:), 'LineWidth', 1)
        plot3(RSave(1,1,i) / R_e, RSave(1,2,i) / R_e, RSave(1,3,i) / R_e,...
            '.', 'color', colors(colorI,:), 'MarkerSize', 10)
        hold on
    end
end

%% DECODE TLE DATA
% ... (plotearth and simStart definitions are fine) ...

% Import the TLE data. 
for k = 1:length(filenames)
    % Get the orbital elements. 
    [coe] = two_line_elem_conv(horzcat(filenames{k}, '.txt'), 'all');
    
    % Find latest epoch (all others can be run up from here).
    coeDateNums = datenum(coe.date);
    [val, ind] = min(coeDateNums);
    
    % Define the max time from simStart.
    a = coe.a(1);
    n = sqrt(mu/a^3);
    tFinal = 2*pi/n*1.1/3600/24; % This gives us just over 1 orbit. 
    
    % Create a time vector. 
    tSim = linspace(simStart, simStart + tFinal, 200);
    
    % Allocate space for plots AND Doppler data
    RSave = NaN(length(tSim), 3, length(coeDateNums));
    DopplerShift_Hz = NaN(length(tSim), length(coeDateNums));
    Elevation_deg = NaN(length(tSim), length(coeDateNums));
    
    % Earth rotation vector (for ECI-to-ECEF velocity conversion)
    omega_E_vec = [0; 0; omega_E];

    % Run through all of the satellites in the TLE file 
    for i = 1:length(coeDateNums) % For each satellite
        
        % Get this satellite's ECI orbital elements
        a = coe.a(i);
        n_sat = sqrt(mu / a^3); % [rad/s]
        e = coe.e(i);
        inc = coe.i(i) * pi / 180;
        RAAN_eci = coe.RAAN(i) * pi / 180;
        omega = coe.omega(i) * pi / 180;
        M_epoch = coe.M(i) * pi / 180;
        t_epoch = coeDateNums(i);
            
        for j = 1:length(tSim) % For each time step
            % --- 1. Propagate Satellite to tSim(j) ---
            dt_sec = (tSim(j) - t_epoch) * 24 * 3600;
            M = M_epoch + n_sat * dt_sec; 
            
            % --- 2. Get Satellite ECI State ---
            % Assumes COE2RV outputs [R_eci_col, V_eci_col]
            % ** You must modify your code to get Velocity **
            [R_sat_eci, V_sat_eci] = COE2RV(a, e, inc, RAAN_eci, omega, M);
            
            % Save ECI position for plotting
            % NOTE: Your original GMST correction was problematic.
            % For a correct plot, you should rotate R_sat_eci by gmst_j
            % *before* plotting, or plot in ECI and rotate the Earth.
            % For simplicity, I'll just save the raw ECI vector.
            RSave(j,:,i) = R_sat_eci'; 

            % --- 3. Get GMST at current time ---
            gmst_j_rad = utc2gmst(datevec(tSim(j)));
            
            % --- 4. Convert Satellite State to ECEF ---
            % ECI to ECEF Rotation Matrix
            R_z = [ cos(gmst_j_rad), sin(gmst_j_rad), 0;
                   -sin(gmst_j_rad), cos(gmst_j_rad), 0;
                   0,                0,               1 ];
            
            R_sat_ecef = R_z * R_sat_eci;
            % Transport theorem: V_ecef = R_z * V_eci - (omega_E x R_ecef)
            V_sat_ecef = R_z * V_sat_eci - cross(omega_E_vec, R_sat_ecef);
            
            % --- 5. Calculate Relative Vectors (in ECEF) ---
            R_rel = R_sat_ecef - R_ue_ecef;
            V_rel = V_sat_ecef - V_ue_ecef;
            
            range = norm(R_rel);
            R_rel_unit = R_rel / range;
            
            % --- 6. Calculate Elevation Angle (Visibility) ---
            % Convert relative ECEF vector to local ENU (East-North-Up)
            R_enu = [ -sin(lon_ue_rad),                cos(lon_ue_rad),               0;
                      -sin(lat_ue_rad)*cos(lon_ue_rad), -sin(lat_ue_rad)*sin(lon_ue_rad), cos(lat_ue_rad);
                       cos(lat_ue_rad)*cos(lon_ue_rad),  cos(lat_ue_rad)*sin(lon_ue_rad), sin(lat_ue_rad) ];
            
            R_rel_enu = R_enu * R_rel;
            
            % El = atan2(Up, sqrt(East^2 + North^2))
            el_rad = atan2(R_rel_enu(3), sqrt(R_rel_enu(1)^2 + R_rel_enu(2)^2));
            Elevation_deg(j, i) = rad2deg(el_rad);
            
            % --- 7. Calculate Doppler (if visible) ---
            if Elevation_deg(j, i) > 5.0  % 5-degree elevation mask
                % Radial velocity (projection of V_rel onto R_rel)
                v_r = dot(V_rel, R_rel_unit);
                
                % Doppler shift
                DopplerShift_Hz(j, i) = f_carrier * (v_r / c_light);
            else
                DopplerShift_Hz(j, i) = NaN; % Not visible
            end
        end
    end
    
    % Plot the orbit (This part of your code is fine)
    % Compute GMST at simStart for the plot
    GMST = utc2gmst(datevec(simStart)); % [rad]
    
    for i = 1:length(coeDateNums)
        % Create the rotation matrix for the plot snapshot time
        R_z_plot = [ cos(GMST), sin(GMST), 0;
                    -sin(GMST), cos(GMST), 0;
                     0,         0,        1 ];
        
        % Rotate all ECI positions to ECEF *at the simStart time*
        R_plot_ecef = (R_z_plot * RSave(:,:,i)')'; 
        
        colorI = k;
        plot3(R_plot_ecef(:,1) / R_e, R_plot_ecef(:,2) / R_e, R_plot_ecef(:,3) / R_e,...
            'color', colors(colorI,:), 'LineWidth', 1)
        plot3(R_plot_ecef(1,1) / R_e, R_plot_ecef(1,2) / R_e, R_plot_ecef(1,3) / R_e,...
            '.', 'color', colors(colorI,:), 'MarkerSize', 10)
        hold on
    end
end
% ... (Figure saving code is fine) ...

%% PLOT DOPPLER INFORMATION
figure;
t_hours = (tSim - tSim(1)) * 24;
plot(t_hours, DopplerShift_Hz / 1e3); % Plot in kHz
grid on;
xlabel('Time from Sim Start (hours)');
ylabel('Doppler Shift (kHz)');
title(sprintf('Doppler Shift for %s at %.2f°N, %.2f°W (f_0=%.2f MHz)', ...
    filenames{1}, lat_ue_deg, -lon_ue_deg, f_carrier/1e6));
legend(coe.satnames, 'Location', 'best');

%% PLOT ELEVATION
figure;
plot(t_hours, Elevation_deg);
grid on;
xlabel('Time from Sim Start (hours)');
ylabel('Elevation Angle (degrees)');
title(sprintf('Elevation for %s at %.2f°N, %.2f°W', ...
    filenames{1}, lat_ue_deg, -lon_ue_deg));
legend(coe.satnames, 'Location', 'best');

%% SAVE FIGURE

% If you want a black background set to 'off', otherwise set to 'on' or
% just comment this out.
%set(gcf, 'InvertHardCopy', 'on');

% Set the view angle of the figure. 
view([-20, 9])

% Reset the zoom. 
% zoom reset
zoom(1.25)

% Turn off axis clipping. 
ax = gca;               
ax.Clipping = 'off';    

% Export the figure. 
exportfig(gcf,horzcat(filenames{1},'Multi.tiff'),'height',6,'width',9,'fontsize',16,'LineWidth',10,'resolution',220);
