function [Fs,humidity,temperature,order,H_length,H_filename,air_F,smooth_F,HPass_F,plot_F2,plot_F3,dist_F,alpha_F...
   ,c,room_size,receiver,sensor,sensor_space,sensor_dir,S_No,source_polar,F_abs,A]=excel_setup(filename)
%Usage: [Fs,humidity,temperature,order,H_length,H_filename,air_F,smooth_F,HPass_F,plot_F2,plot_F3,dist_F,alpha_F...
%    ,c,room_size,receiver,sensor,sensor_space,sensor_dir,S_No,source_polar,F_abs,A]=excel_setup(filename);
% EXCEL_SETUP.M reads the Excel spreadsheet file SETUP.XLS to obtain values for setting up roomsim.
% It calls GetExcelSetup.m, a version of the standard Matlab XLSREAD.M function modified to handle only text data
%   and to allow compilation to stand-alone code.
%---------------------------------------------------------------------------- 
% Copyright (C) 2003  Douglas R Campbell
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program in the MATLAB file roomsim_licence.m ; if not,
%  write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
%---------------------------------------------------------------------------------
% Functions called: GetExcelSetup.m, 

global hrir_l hrir_r; % Globals for CIPIC impulse responses left and right
global H3D; % Global for MIT impulse response data

sensor_dir=[]; % Declare empty array
error_title='Excel_setup: Error';

%------------- Get the data from the setup Excel spreadsheet -----------------------
[Param_text] = GetExcelSetup(filename,'single values');
[sources_text] = GetExcelSetup(filename,'sources');
[dir_text] = GetExcelSetup(filename,'sensor dir');
[absorption_text] = GetExcelSetup(filename,'surface absorption');

%----------------------- Simulation control parameters ------------------------------
Fs=str2num(Param_text{2,2}); % Sampling frequency (Hz)
humidity=str2num(Param_text{3,2}); % Relative humidity of air (%) (Used to calculate air absorption coefficient "m", valid range 20%< h <70%)
temperature=str2num(Param_text{4,2}); % Temperature of air (deg C) (Used to calculate speed of sound (m/s))
order=str2num(Param_text{5,2}); % If -ve then value computed in make_Roomsim is used, else value supplied here is used (limits order of reflections computed)
H_length=str2num(Param_text{6,2}); % If -ve then H_length is later set = RT60, else value supplied here is used.

%--------------------------------- Output destination --------------------------------
H_filename=Param_text{7,2}; %Output filename for impulse response.

%--------------------------------------- Flags ---------------------------------------
air_F=logical(str2num(Param_text{8,2})); % false = no absorption due to air, true = air absorption is present.
smooth_F=logical(str2num(Param_text{9,2})); %false = no smoothing filter applied, true = smoothing filter used.
HPass_F=logical(str2num(Param_text{10,2})); %false = no High-Pass filter applied, true = 100 Hz High-Pass filter applied.
plot_F2=logical(str2num(Param_text{11,2})); % false = no plot, true = 2D-plan, shows image rooms on xy plane.
plot_F3=logical(str2num(Param_text{12,2})); % false = no plot, true = 3D-plot, rotatable.

%--------- Room dimensionss --------------
Lx=str2num(Param_text{13,2}); % Height
Ly=str2num(Param_text{14,2}); % Width
Lz=str2num(Param_text{15,2}); % Length

% Receiver reference point [xp,yp,zp](m), if Head, it is mid-point of the inter_aural axis
xp=str2num(Param_text{16,2}); % x coordinate of receiver refrence point
yp=str2num(Param_text{17,2}); % y coordinate of receiver refrence point
zp=str2num(Param_text{18,2});% z e.g. Typical height above floor of ears of seated human subject = 1.2 m

% Sensor details
sensor=Param_text{19,2};
sensor_space=str2num(Param_text{20,2}); % Sensor separation (if head it is implicit in the HRIR data)

% File control values
dir_ch=Param_text{21,2};
ext=Param_text{22,2};

%-------------- Identify path to filename for extraction of hrir from MIT Kemar data base -------------
MIT_root=Param_text{23,2}; % Root directory of the MIT Kemar data, an immediate sub-directory of the Roomsim directory
MIT_subdir=Param_text{24,2};
MIT_filename=Param_text{25,2};

%-------------- Identify path to filename for extraction of hrir from CIPIC data base -----------------
S_No=Param_text{26,2}; % CIPIC subject number, format '&&&' (e.g. '021' Kemar with small pinnae)
CIPIC_root=Param_text{27,2}; % Root directory of the CIPIC HRTF data, an immediate sub-directory of the Roomsim directory
CIPIC_subdir=Param_text{28,2};
CIPIC_filename=Param_text{29,2};

%------------- Debug Flags --------------
dist_F=logical(str2num(Param_text{30,2})); % false = no distance effect, true = 1/R attenuation with distance applied.
alpha_F=logical(str2num(Param_text{31,2})); % false = fixed transparent surfaces for Room Geometry plot, true = (surface opacity = reflectivity)

% Sound source position(s) relative to receiver reference point (Head or sensor(s)).
limit=size(sources_text,2); %Get the number of columns in the sheet (allows user to add sources)
for k=1:limit-1, %
    if ~(isempty(sources_text{2,k+1})||isempty(sources_text{3,k+1})||isempty(sources_text{4,k+1})), % Prevent reading an empty cell
        R_s(1,k)=str2num(sources_text{2,k+1}); % Row vector of radial distance(s) of source(s) from head (m)
        alpha(1,k)=str2num(sources_text{3,k+1}); % Row vector of azimuth(s) of sources -180< alpha < 180 (deg) NB +ve is ACW on xy plane
        beta(1,k)=str2num(sources_text{4,k+1}); % Row vector of elevation(s) of sources -90< beta < 90 (deg).
    end;
end;

switch sensor % Selection of sensor type 
    case {'one_mic'}
        % Set up sensor directionality  (Omni is -180,180,-90,90 deg)
        min_azim_sensor=str2num(dir_text{2,2}); % Minimum azimuth seen by sensor
        max_azim_sensor=str2num(dir_text{3,2}); % Maximum azimuth seen by sensor
        min_elev_sensor=str2num(dir_text{4,2}); % Minimum elevation seen by sensor
        max_elev_sensor=str2num(dir_text{5,2}); % Maximum elevation seen by sensor
        sensor_dir=[min_azim_sensor; max_azim_sensor; min_elev_sensor; max_elev_sensor]; % Pack up sensor directionality one column per sensor

    case {'two_mic'}
        % Set up sensor directionality [Left Right] (Omni is [-180 -180],[180 180],[-90 -90],[90 90] deg)
        min_azim_sensor=[str2num(dir_text{6,2}) str2num(dir_text{6,3})];
        max_azim_sensor=[str2num(dir_text{7,2}) str2num(dir_text{7,3})];
        min_elev_sensor=[str2num(dir_text{8,2}) str2num(dir_text{8,3})];
        max_elev_sensor=[str2num(dir_text{9,2}) str2num(dir_text{9,3})];
        sensor_dir=[min_azim_sensor; max_azim_sensor; min_elev_sensor; max_elev_sensor]; % Pack up sensor directionality one column per sensor

    case 'mithrir' % Extract hrir from MIT Kemar data base
        if (exist(MIT_root,'dir')==7) % Check path name to the MIT Kemar directory (folder)
            MIT_file = [MIT_root dir_ch MIT_subdir dir_ch MIT_filename ext]; %Form pathname to MIT kemar file
            load (MIT_file); % Load the MIT Kemar file into the workspace (creates H3D)
        else,
            h=errordlg('Directory (folder) MIT Kemar data file (hrir_final.mat) not found, exiting from roomsim','error_title');
            beep;
            uiwait(h);
            return;
        end;
        
    case 'cipicir' % Extract hrir from CIPIC data base
        if (exist(CIPIC_root,'dir')==7),
            CIPIC_subpath=[CIPIC_subdir '_' S_No];
            CIPIC_file = [CIPIC_root dir_ch CIPIC_subpath dir_ch CIPIC_filename ext]; %Form pathname to CIPIC subject file
            load(CIPIC_file); % Load the subject file into the workspace (creates hrir_l and hrir_r)
        else,
            h=errordlg('Directory (folder) CIPIC not found, exiting from roomsim','error_title');
            beep;
            uiwait(h);
            return;
        end;
end

%---------------- Set the room surface absorptions ----------------------
for k=1:6,
    F_abs(k)=str2num(absorption_text{2,k+1});
    Ax1(k)=str2num(absorption_text{3,k+1}); % Absorption of wall in x=0 plane (behind Kemar in plan)
    Ax2(k)=str2num(absorption_text{4,k+1}); % Absorption of wall in x=Lx plane (front in plan)
    Ay1(k)=str2num(absorption_text{5,k+1}); % Absorption of wall in y=0 plane (right in plan)
    Ay2(k)=str2num(absorption_text{6,k+1}); % Absorption of wall in y=Ly plane (left in plan)
    Az1(k)=str2num(absorption_text{7,k+1}); % Absorption of floor i.e. z=0 plane
    Az2(k)=str2num(absorption_text{8,k+1}); % Absorption of ceiling i.e. z=Lz plane
end;
%---------------- End of surface absorption set up ---------------

c=round(331*sqrt(1+0.0036*temperature)); % Calculate speed of sound (m/s) as function of temperature

% Pack various parameters for the room simulation
room_size=[Lx;Ly;Lz]; % Pack up room dimensions into column vector.
receiver=[xp;yp;zp]; % Pack up receiver (listener's head) coordinates into column vector.
source_polar=[R_s;alpha;beta]; % Pack up source(s) coordinates into column vector.
A=[Ax1' Ax2' Ay1' Ay2' Az1' Az2']; % Pack up column vectors (NB Transposition) of absorption coefficients in array A
%----------------------------------------- End of excel_setup.m ----------------------------------------
