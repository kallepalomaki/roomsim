function [pulse, elerr, azerr] = getNearestMITpulse(elev,azim,H3D);
%Usage: [pulse, elerr, azerr] = getNearestMITpulse(elev,azim,H3D);
%Retrieves from H3D, the impulse response (a cell array) that is closest to the specified elevation and azimuth (deg)
%----------------------------------------------------------------------------- 
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
%-----------------------------------------------------------------------------------
%Functions called: 

if nargin < 1
    fprintf('Format: [pulse, azerr, elerr] = getNearestMITpulse(azim,elev,H3D)\n');
    return;
end

% quan.m Quantises the elevations and azimuths to nearest angle in the Kemar set
[elevation,azimuth] = quan(elev,-azim); % NB azimuth for Kemar measurements was +ve CW so invert sign.

elevations = [-40 -30 -20 -10  0 10 20 30 40 50 60 70 80 90;
	           56  60  72  72 72 72 72 60 56 45 36 24 12  1];
  
%Quantize elevation
[elerr el]=min(abs(elevations(1,:)-elevation)); % Find the index el of the value within elevations which is the minimum distance from ele
elevation=elevations(1,el);

%Quantize azimuth
if azimuth < 0 
  azimuth=360+azimuth;
end

azim_incr = 360/elevations(2,el);
azimuths = round([0:azim_incr:360-azim_incr]);
[azerr az]=min(abs(azimuths-azimuth)); % Find the index az of the value within azimuths which is the minimum distance from azimuth
azimuth=azimuths(az);

if azimuth>180
  azimuth=azimuth-360; % Restrict to range -180 < azimuth <= 180
end

pulse = H3D{el,az}; % Extract the required hrir L&R pair (Left is first column ie pulse(:,1) )
%--------------------- End of getNearestMITpulse.m -----------------------------
