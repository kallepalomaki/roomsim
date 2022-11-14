function schematic_head(coordinates,head_width,skin)
% Usage: schematic_head(coordinates,head_width,skin);  Draw a schematic head with eyes, 
% nose and ears, head_width, ear separation for head (m), skin is the colormap for sphere.
%------------------------------------------------------------------------------------ 
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
%--------------------------------------------------------------------------------------
%Functions called: c_tetrahedron3.m

% Unpack the reference coordinates locating the head position
x=coordinates(1);
y=coordinates(2);
z=coordinates(3);

scale_F=0.5*head_width;
[XX YY ZZ]=sphere; % Data for drawing head
cm=colormap(skin); % Skin colour
cmin=min(min(z+scale_F*ZZ));
cmax=max(max(z+scale_F*ZZ));
caxis([cmin cmax]); % Scale the colormap
C_Faces=0.65*(cmin+cmax); %Colour of nose
C_Faces3=[1 0 0;0 0 1;0 0 0]; % Nose colours red (R), blue (L), black (base)
C_Edge='k'; % Edge colour of tetrahedral nose
hold on;
surf(x+scale_F*XX, y+scale_F*YY, z+1.25*scale_F*ZZ,'EdgeColor','none');
%Draw the ears
surf(x+0.15*scale_F*XX, y+0.98*scale_F+0.15*scale_F*YY, z+0.15*scale_F*ZZ,'EdgeColor','b');
surf(x+0.15*scale_F*XX, y-0.98*scale_F+0.15*scale_F*YY, z+0.15*scale_F*ZZ,'EdgeColor','r');
% Draw the nose
c_tetrahedron3(x+0.9*scale_F,y,z+0.1*scale_F,0.15*scale_F,C_Faces3,1,C_Edge);
% Draw the eyes
surf(x+0.7*scale_F+0.15*scale_F*XX, y+0.5*scale_F+0.15*scale_F*YY, z+0.2*scale_F+0.15*scale_F*ZZ,'EdgeColor','w'); % Draw sclera
surf(x+0.7*scale_F+0.15*scale_F*XX, y-0.5*scale_F+0.15*scale_F*YY, z+0.2*scale_F+0.15*scale_F*ZZ,'EdgeColor','w');
surf(x+0.8*scale_F+0.05*scale_F*XX, y+0.55*scale_F+0.05*scale_F*YY, z+0.2*scale_F+0.05*scale_F*ZZ,'EdgeColor','b'); % Draw pupil
surf(x+0.8*scale_F+0.05*scale_F*XX, y-0.55*scale_F+0.05*scale_F*YY, z+0.2*scale_F+0.05*scale_F*ZZ,'EdgeColor','r');
axis equal;
%--------------------------------- End of schematic_head.m -------------------------------
