function [pulse, azerr, elerr] = getNearestUCDpulse(azimuth,elevation,h3D);
%Usage: [pulse, azerr, elerr] = getNearestUCDpulse(azimuth,elevation,h3D);
%Retrieves the impulse response from h3D that is closest to the
%specified azimuth and elevation (in degrees0
%-----------------------------------------------------------------------
% Acknowledgement: This function is only marginally changed from that
% published by CIPIC in their "hrir_data_documentation.pdf"
%----------------------------------------------------------------------
% Functions called: pvaldeg.m

if nargin < 1
    fprintf('Format: [pulse, azerr,elerr] = getNearestUCDpulse(azimuth,elevation,h3D)\n');
    return;
end

azimuth = pvaldeg(azimuth);
if (azimuth < -90) | (azimuth > 90)
    error('Invalid azimuth');
end
elevation = pvaldeg(elevation);

elmax = 50;
elindices = 1:elmax;
elevations = -45 + 5.625*(elindices-1);
el = round((elevation+45)/5.625 + 1);
el = max(el,1);
el = min(el,elmax); %Index for required elevation
elerr = pvaldeg(elevation - elevations(el)); %Error in using this indexed value

azimuths = [-80 -65 -55 -45:5:45 55 65 80];
[azerr, az] = min(abs(pvaldeg(abs(azimuths - azimuth))));%Error and Index for required azimuth

pulse = squeeze(h3D(az,el,:));
%--------------- End of getNearestUCDpulse.m --------------------------------------------
