function [RT60, MFP]= Norris_Eyring(c,room_size,A,F_abs,m_air)
%Usage: [RT60, MFP]= Norris_Eyring(c,room_size,A,F_abs,m_air);
% Calculate a room reverberation time (RT60) using the Norris-Eyring estimate
%Ref. http://www.teicontrols.com/notes/AcousticsEE363N/EngineeringAcoustics.pdf
%-------------------------------------------------------------------------------------- 
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
%------------------------------------------------------------------------------------
%Functions called: 

Lx=room_size(1);
Ly=room_size(2);
Lz=room_size(3);
V_room = Lx*Ly*Lz; % Volume of room m^3
Sxz=Lx*Lz;
Syz=Ly*Lz;
Sxy=Lx*Ly;
S =2*(Sxz + Syz + Sxy); % Total area of shoebox room surfaces
Se = Syz*(A(:,1) + A(:,2)) + Sxz*(A(:,3) + A(:,4)) + Sxy*(A(:,5) + A(:,6)); %Total absorption of room surfaces at each frequency
a_bar = Se/S; % Mean absorption of room surfaces
m = mean(m_air); % Mean absorption of air averaged across frequency.
MFP = 4*V_room/S; % Mean Free Path (Average distance between succesive reflections) (Ref A4)

% Norris-Eyring estimate adjusted for air absorption
if abs(1-a_bar) < eps % Detect anechoic case and force RT60 all zeros.
    RT60 = zeros(size(F_abs));
else %Norris-Eyring reverberation time estimate
    RT60 = (55.25/c)*V_room./(4*m_air'*V_room-S*log(1-a_bar));
end
%-------------------------- End of Norris-Eyring.m ---------------------------------------------
