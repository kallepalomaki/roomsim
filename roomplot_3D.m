function roomplot_3D(room_size,source_xyz,receiver,sensor,sensor_xyz,sensor_dir,source,p_isource,fig_loc)
%Usage: roomplot_3D(room_size,source_xyz,receiver,sensor,sensor_xyz,sensor_dir,source,p_isource);
% Draw a 3D view of the room and its image sources
%--------------------------------------------------------------------------------------------- 
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
%---------------------------------------------------------------------------------------------------
% Functions called: 

%*** Global Declarations ***
global curr_fig; %Keep track of current figure number

%---------------------------------- Initialisation -----------------------------
invis_h=[]; imagesM_h=[]; imagesL_h=[]; imagesR_h=[]; % Declare arrays as empty to avoid error if no such to be plotted
[L_colormap, R_colormap, LR_colormap]=isource_colormaps; % Get the user defined colormaps L, R & LR.
n_sources=size(p_isource,2); %Number of primary sources
n_images=size(p_isource,1); % Number of image sources

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

button=[];
while isempty(button) % Loop until valid choice has been made (disable close window button)
    button = questdlg('Show "invisible" image sources?','Roomplot_3D','Yes','No','No');
    beep;
end; % Of while loop to disable inappropriate close window button operation
if strcmp(button,'Yes')
    vis_F=true; % Show invisibles/inaudibles
else
    vis_F=false; % Hide invisibles/inaudibles
end;

%------------- Specific setups for each sensor case -------------------------------
switch sensor % ---------- Identify receiver system---------------------
    case {'one_mic','two_mic'}
        % Unpack sensor directionality
        min_azim_sensor=sensor_dir(1,:);
        max_azim_sensor=sensor_dir(2,:);
        min_elev_sensor=sensor_dir(3,:);
        max_elev_sensor=sensor_dir(4,:);
end;

%----------------------- Display a 3D plot of the room, receiver and the image sources ---------------
curr_fig=get(0,'CurrentFigure'); %Get handle of current figure
if isempty(curr_fig),
    figure(1); %Start first figure
else
    figure(curr_fig+1); %Start new figure
end
set(gcf,'position',fig_loc); % Put figure at this location
clf; % Make sure it's a clean sheet
back_col=[0.7 0.7 0.7]; % Background colour light grey
whitebg(back_col); % Set figure background colour
hold on

%------------------------------- Draw the room in 3D ------------------------
alph=[62; 20; 20; 20; 20; 20]; %Transparency of faces fixed NB order here is Az1,Az2,Ax1,Ax2,Ay1,Ay2
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

%Plot the location of the image sources
x_size=[]; y_size=[]; z_size=[]; % Clear accumulators for 3Dplot sizing
switch sensor % ---------- Identify receiver system---------------------
    case {'one_mic'}
        sensorM_h=plot3(sensor_xyz(1),sensor_xyz(2),sensor_xyz(3),'o','MarkerEdgeColor','k','MarkerFaceColor','m'); % Position of sensor       
        cm_LR = LR_colormap; % Get the 60 colormap RGB values
        for ps=1:n_sources % For each parent source
            for is=1:n_images % For each image source (with delays <H_length or impulse response peaks >-60dB)
                if p_isource(is,ps)>0 % Plot only "visible" sources
                    colour_index=round(60+20*log10(p_isource(is,ps,1))); % Form Index from bottom (1) to top (60) of the colormap
                    if colour_index <1, colour_index=1; elseif colour_index >60, colour_index=60; end;
                    mark_colour= cm_LR(colour_index,:);
                    imagesM_h=plot3(source(1,is,ps),source(2,is,ps),source(3,is,ps),'+','MarkerEdgeColor',mark_colour); %Plot the visible image source positions (xyz,ImageNo,Source)
                elseif vis_F % Plot "invisible" image sources in black
                    invis_h=plot3(source(1,is,ps),source(2,is,ps),source(3,is,ps),'k.'); %Plot positions for invisible image source(xyz,ImageNo,ParentSource)
                end;    
                temp_x= source(1,is,ps); temp_y= source(2,is,ps); temp_z= source(3,is,ps);
                x_size=[x_size temp_x]; y_size=[y_size temp_y]; z_size=[z_size temp_z]; % Accumulate size for plotting
            end; % For each image source
        end; % For each parent source
        
    case {'two_mic','mithrir','cipicir'} % Assumes two-sensor inter-sensor axis is parallel with y axis and lies in z=constant plane 
        plot3(sensor_xyz(1,:),sensor_xyz(2,:),sensor_xyz(3,:),'k',receiver(1),receiver(2),receiver(3),'.k'); %%Draw line between L & R sensors and mark the receiver reference point
        sensorL_h=plot3(sensor_xyz(1,1),sensor_xyz(2,1),sensor_xyz(3,1),'o','MarkerEdgeColor','k','MarkerFaceColor','b'); % Position of L sensor     
        sensorR_h=plot3(sensor_xyz(1,2),sensor_xyz(2,2),sensor_xyz(3,2),'o','MarkerEdgeColor','k','MarkerFaceColor','r'); % Position of R sensor       

        cm_L = L_colormap; % Get the 60 colormap RGB values for left sensor only
        cm_R = R_colormap; % Get the 60 colormap RGB values for right sensor only
        for ps=1:n_sources % For each parent source
            for is=1:n_images % For each image source
                if p_isource(is,ps,1)>0 % Plot sources "visible" to L sensor
                    colour_index=round(60+20*log10(p_isource(is,ps,1))); %Index from bottom (1) to top (60) of the colormap
                    if colour_index <1, colour_index=1; elseif colour_index >60, colour_index=60; end;
                    mark_colour= cm_L(colour_index,:);
                    imagesL_h=plot3(source(1,is,ps),source(2,is,ps),source(3,is,ps),'+','MarkerEdgeColor',mark_colour); %Plot the visible image source positions (+) (xy,ImageNo,Source)
                end;
                if p_isource(is,ps,2)>0 % Plot sources "visible" to R sensor
                    colour_index=round(60+20*log10(p_isource(is,ps,2))); %Index from bottom (1) to top (60) of the colormap
                    if colour_index <1, colour_index=1; elseif colour_index >60, colour_index=60; end;
                    mark_colour= cm_R(colour_index,:);
                    imagesR_h=plot3(source(1,is,ps),source(2,is,ps),source(3,is,ps),'s','MarkerEdgeColor',mark_colour); %Plot the visible image source positions (square) (xy,ImageNo,Source)
                end;
                if vis_F && (p_isource(is,ps,1)== 0 && p_isource(is,ps,2) == 0) % Plot "invisible" image sources in black
                    invis_h=plot3(source(1,is,ps),source(2,is,ps),source(3,is,ps),'k.'); %Plot positions of invisible image sources (xyz,ImageNo,ParentSource)
                end
                temp_x= source(1,is,ps); temp_y= source(2,is,ps); temp_z= source(3,is,ps);
                x_size=[x_size temp_x]; y_size=[y_size temp_y]; z_size=[z_size temp_z]; % Accumulate size for plotting
            end; % For each image source
        end; % For each parent source
        
    otherwise
        disp('Unknown receiver set up at freq response plot');
        return
end; % Switch sensor

%Scale the plot to show all image sources and the complete room.
% xmin= min(min(source(1,:,:))); xmax= max(max(source(1,:,:)));
% ymin= min(min(source(2,:,:))); ymax= max(max(source(2,:,:)));
% zmin= min(min(source(3,:,:))); zmax= max(max(source(3,:,:)));
xmin= min([x_size -Lx])-Lx; xmax= max([x_size Lx]+Lx);
ymin= min([y_size -Ly])-Ly; ymax= max([y_size Ly]+Ly);
zmin= min([z_size -Lz])-Lz; zmax= max([z_size Lz]+Lz);

axis([xmin xmax ymin ymax zmin zmax]);
xlabel('Length (x) m'); ylabel('Width (y) m'); zlabel('Height (z) m');
title('3D Plot of Room and Image Sources');
axis tight equal ;
view(3); %Set the default 3D view aspect
grid;
switch sensor % ---------- Identify the receiver system and display the legend ---------------------
    case {'one_mic'}
        if vis_F, %Show legend that features "invisible" image sources
            legend([outline_h; source_h; sensorM_h; imagesM_h; invis_h],{'Ay1';'Ax1';'Ay2';'Ax2'...
                    ;'Source(s)';'Receiver';'Images Mono';'"Inaudible"'},-1);
        else
            legend([outline_h; source_h; sensorM_h; imagesM_h],{'Ay1';'Ax1';'Ay2';'Ax2'...
                    ;'Source(s)';'Receiver';'Images Mono'},-1);
        end;
    case {'two_mic','mithrir','cipicir'} % Assumes two-sensor inter-sensor axis is parallel with y axis and lies in z=constant plane
        if vis_F, %Show legend that features "invisible" image sources
            legend([outline_h; source_h; sensorL_h; sensorR_h; imagesL_h; imagesR_h; invis_h],{'Ay1';'Ax1';'Ay2';'Ax2'...
                    ;'Source(s)';'Receiver L';'Receiver R';'Images L';'Images R';'"Inaudible"'},-1);
        else
            legend([outline_h; source_h; sensorL_h; sensorR_h; imagesL_h; imagesR_h],{'Ay1';'Ax1';'Ay2';'Ax2'...
                    ;'Source(s)';'Receiver L';'Receiver R';'Images L';'Images R'},-1);
        end;
end;

hold off;
drawnow; %Force completion of all previous figure drawing before continuing
%------------------------------ End of roomplot_3D.m ----------------------------------------------
