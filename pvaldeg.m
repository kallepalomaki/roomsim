function pangle = pvaldeg(theta)
%Usage: pangle = pvaldeg(theta); Maps angle theta (in degrees) into the range (-90, 270)
%----------------------------------------------------------------------
% Ackowledgement: This function is only marginally changed from that
% published by CIPIC in their "hrir_data_documentation.pdf"
%-----------------------------------------------------------------------
% Functions called: 

if nargin < 1
    fprintf('Format: pangle = pvaldeg(theta)\n');
    return
end
rtd=180/pi; %radians to degrees factor
theta_rad=theta*pi/180;
pangle = atan2(sin(theta_rad),cos(theta_rad))*rtd;

if pangle < -90
    pangle = pangle + 360;
end
%---------------------- End of pvaldeg.m -------------------------