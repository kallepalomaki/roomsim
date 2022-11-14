function roomplot(room_size,source_xyz,receiver,sensor,sensor_xyz,head_width,alph,fig_loc)
%Usage: roomplot(room_size,source_xyz,receiver,sensor,sensor_xyz,head_width,alph);
%Display the room geometry with receiver and source(s) locations as a 3D plot for visual confirmation
%-------------------------------------------------------------
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
%---------------------------------------------------------------------
% Functions called: schematic_head.m

%*** Global Declarations ***
global curr_fig; %Keep track of current figure number

% Unpack room dimensions
Lx=room_size(1); % Length
Ly=room_size(2); % Width
Lz=room_size(3); % Height

%Unpack source(s) coordinates
x=source_xyz(1,:);
y=source_xyz(2,:);
z=source_xyz(3,:);

%Unpack receiver reference coordinates (Head or sensor(s))
xp=receiver(1);
yp=receiver(2);
zp=receiver(3);

curr_fig=get(0,'CurrentFigure'); %Get handle of current figure
if isempty(curr_fig),
    figure(1); %Start first figure
else
    figure(curr_fig+1); %Start new figure
end
set(gcf,'position',fig_loc); % Put figure at this location
clf ; %make sure it's a clean sheet
back_col=[0.7 0.7 0.7]; % Background colour light grey
whitebg(back_col); % Set figure background colour
hold on;

%--------------- Draw the room in 3D ------------------------
% Rearrange face transparency (alpha) order from Ax1,Ax2,Ay1,Ay2,Az1,Az2 to Az1,Az2,Ax1,Ax2,Ay1,Ay2
temp1=alph(5); % Save Az1
temp2=alph(6); % Save Az2
alph(3:6)=alph(1:4); % Copy Ax1,Ax2,Ay1,Ay2 to bottom four locations
alph(1)=temp1; % Put Az1 at top
alph(2)=temp2; % Put Az2 at second top

V1=[0 0 0]; V2=[Lx 0 0]; V3=[0 Ly 0]; V4=[0 0 Lz]; V5=[Lx Ly 0]; V6=[Lx 0 Lz]; V7=[0 Ly Lz]; V8=[Lx Ly Lz]; %Room verteces [x y z]
vertex_list=[V1;V2;V3;V4;V5;V6;V7;V8];
faces=[V1 V2 V5 V3; V4 V6 V8 V7; V1 V4 V7 V3; V2 V6 V8 V5; V1 V2 V6 V4; V3 V5 V8 V7]; % Az1,Az2,Ax1,Ax2,Ay1,Ay2
vertex_connection=[1 2 5 3;4 6 8 7;1 4 7 3;2 6 8 5;1 2 6 4;3 5 8 7]; % Order of connection of vertices
CData=[0.5 0.5 0.5;1 1 1;1 0 1;0 1 0;1 0 0;0 0 1]; % Face colours grey,white,magenta,green,red,blue
patch('Faces',vertex_connection,'Vertices',vertex_list,'AlphaDataMapping','direct','FaceAlpha','flat','FaceVertexAlphaData',alph,'FaceColor','flat','FaceVertexCdata',CData,'EdgeColor','none');    
% Outline the room edges
outline_h=plot3([0 Lx],[0 0],[0 0],'r',[0 0],[0 Ly],[0 0],'m',[0 Lx],[Ly Ly],[0 0],'b',[Lx Lx],[0 Ly],[0 0],'g'); % Outline the floor (Solo for Legend)
plot3([Lx Lx],[0 0],[0 Lz],'r',[0 0],[0 0],[0 Lz],'r',[0 0],[Ly Ly],[0 Lz],'b',[Lx Lx],[Ly Ly],[0 Lz],'b'...%  Outline the vertical edges
    ,[0 Lx],[0 0],[Lz Lz],'r',[0 0],[0 Ly],[Lz Lz],'m',[0 Lx],[Ly Ly],[Lz Lz],'b',[Lx Lx],[0 Ly],[Lz Lz],'g'); % Outline the ceiling 
%---------------------------------------- End of draw room -------------------------------------------

source_h=plot3(x,y,z,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); % Position of source(s)
switch sensor % ---------- Identify and position the receiver system---------------------
    case {'one_mic'}
        sensorM_h=plot3(sensor_xyz(1),sensor_xyz(2),sensor_xyz(3),'o','MarkerEdgeColor','k','MarkerFaceColor','m'); % Position of sensor       
    case {'two_mic'} % Assumes two-sensor inter-sensor axis is parallel with y axis and lies in z=constant plane
        plot3(sensor_xyz(1,:),sensor_xyz(2,:),sensor_xyz(3,:),'k'); %Draw line between L & R sensors
        plot3(receiver(1),receiver(2),receiver(3),'.k'); % Mark the receiver reference point
        sensorL_h=plot3(sensor_xyz(1,1),sensor_xyz(2,1),sensor_xyz(3,1),'o','MarkerEdgeColor','k','MarkerFaceColor','b'); % Position of L sensor     
        sensorR_h=plot3(sensor_xyz(1,2),sensor_xyz(2,2),sensor_xyz(3,2),'o','MarkerEdgeColor','k','MarkerFaceColor','r'); % Position of R sensor       
    case {'mithrir','cipicir'} % Assumes inter-ear axis is parallel with y axis and lies in z=constant plane
        %Draw a schematic head for the Receiver 
        [skin]=skin_colormap; % Schematic head colormap
        schematic_head(receiver,head_width,skin);
        sensorL_h=plot3(sensor_xyz(1,1),sensor_xyz(2,1),sensor_xyz(3,1),'o','MarkerEdgeColor','k','MarkerFaceColor','b'); % Position of L ear     
        sensorR_h=plot3(sensor_xyz(1,2),sensor_xyz(2,2),sensor_xyz(3,2),'o','MarkerEdgeColor','k','MarkerFaceColor','r'); % Position of R ear       
end;
axis tight equal;
xlabel('Length (x) '); ylabel('Width (y) '); zlabel('Height (z) ');
title('Room, Receiver & Source(s) Geometry');
view(3); %Set the default 3D view aspect

switch sensor % ---------- Identify the receiver system and display the legend ---------------------
    case {'one_mic'}
        legend([outline_h; source_h; sensorM_h],{'Ay1';'Ax1';'Ay2';'Ax2';'Source(s)';'Receiver'},-1);
    case {'two_mic','mithrir','cipicir'} % Assumes two-sensor inter-sensor axis is parallel with y axis and lies in z=constant plane
        legend([outline_h; source_h; sensorL_h; sensorR_h],{'Ay1';'Ax1';'Ay2';'Ax2'...
                ;'Source(s)';'Receiver L';'Receiver R'},-1);
end;
hold off
drawnow; %Force completion of all previous figure drawing before continuing
%----------------- End of roomplot.m ---------------------------------
