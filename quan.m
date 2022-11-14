function [elevation,azimuth]=quan(ele,azi)
%Usage: [elevation,azimuth]=quan(ele,azi);
% Quantises the elevation and azimuth for the MIT Kemar HRIR data set
%---------------------------------------------------------------------------- 
% Copyright (C) 2003  Douglas R Campbell and Kalle Palomaki
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
%-------------------------------------------------------------------------------
% Functions called: 

elevs = [-40 -30 -20 -10  0 10 20 30 40 50 60 70 80 90;
	      56  60  72  72 72 72 72 60 56 45 36 24 12  1];
  
%Quantize elevation
[mi e_index]=min(abs(elevs(1,:)-ele)); % Find the value within elevs which is the minimum distance from ele
elevation=elevs(1,e_index);

%Quantize azimuth
if azi < 0 
  azi=360+azi;
end

azim_incr = 360/elevs(2,e_index);
azims = round([0:azim_incr:360-azim_incr]);
[mi a_index]=min(abs(azims-azi)); % Find the value within azims which is the minimum distance from azi
azimuth=azims(a_index);

if azimuth>180
  azimuth=azimuth-360; % Restrict to range -180 < azimuth <= 180
end
%---------------- End of quan.m --------------------------