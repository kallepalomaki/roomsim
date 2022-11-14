function roomplot_2D(c,Fs,room_size,source_xyz,receiver,sensor,sensor_xyz,sensor_dir,H_length,source,p_isource,fig_loc)
%Usage: roomplot_2D(c,Fs,room_size,source_xyz,receiver,sensor,sensor_xyz,sensor_dir,H_length,source,p_isource);
% Display a 2D plot of Room, Image Sources and Image Room Boundaries at height slice_z.
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
%--------------------------------------------------------------------------------------------------
%Functions called: 

%*** Global Declarations ***
global curr_fig; %Keep track of current figure number

%--------------------------- Initialisation ----------------------------
invis_h=[]; imagesM_h=[]; imagesL_h=[]; imagesR_h=[]; % Declare arrays as empty to avoid error if no such to be plotted
[L_colormap, R_colormap, LR_colormap]=isource_colormaps; % Get the user defined colormaps L, R & LR.
n_sources=size(p_isource,2); %Number of primary sources
n_images=size(p_isource,1); % Number of image sources
T=1/Fs;

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
%-------------------------------------------------------------------------

%Dialogue for setting 2D plot slice height
answer={};
while isempty(answer),
    banner = 'Roomplot_2D:';
    prompt = {'Enter height for Slice (Default is source1 height z (m) :'};
    lines = 1;
    def = {num2str(z(1))}; %Default value
    answer = inputdlg(prompt,banner,lines,def);
    beep;
end; % Of while loop to disable inappropriate CANCEL button operation
slice_z=str2num(answer{1});% Height value in slice_z.

%Display the visualisation menu
menu_title='Roomplot_2D: Choose a visualisation';
B1_text='Intensity shown by stem height (Pseudo 3D)';
B2_text='Intensity shown by marker colour (2D Plan)';
beep;
M_VIScase=0;
while M_VIScase==0, % Loop until valid choice has been made (disable close window button)
M_VIScase = menu(menu_title,B1_text,B2_text);
switch M_VIScase
    case 1
        stem_F=1; % Display a Pseudo 3D stem plot
    case 2
        stem_F=0; % Display a 2D colour coded plan view
end;
end; % Of while loop to disable inappropriate close window button operation

button=[];
while isempty(button) % Loop until valid choice has been made (disable close window button)
    button = questdlg('Show "invisible" image sources?','Roomplot_2D','Yes','No','No');
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

%----------------------- 2D Plot of the room and the image rooms---------------
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

%-------------------- Draw the room in 2D plan view ------------------------------
outline_h=plot([0 Lx],[0 0],'r',[0 0],[0 Ly],'m',[0 Lx],[Ly Ly],'b',[Lx Lx],[0 Ly],'g'); % Outline the room

source_h=plot(x,y,'o','MarkerEdgeColor','k','MarkerFaceColor','g'); %Position of source(s)

%Plot the location of the image sources
x_size=[]; y_size=[]; % Clear accumulators for 2Dplot sizing

% Display a 2D intensity or 3D stem plot
switch sensor % ---------- Identify and position the receiver system---------------------
    case {'one_mic'}
        sensorM_h=plot(sensor_xyz(1),sensor_xyz(2),'o','MarkerEdgeColor','k','MarkerFaceColor','m'); % Position of mono sensor       
        cm_LR = LR_colormap; % Get the 60 colormap RGB values
        for ps=1:n_sources, % For each parent source
            for is=1:n_images, % For each image source
                if ( (source(3,is,ps) > slice_z*(0.99) & source(3,is,ps) < slice_z*(1.01)) ), % Display a slice at the chosen height +/- 1%
                    temp_x=[]; temp_y=[];
                    if p_isource(is,ps,1) > 0 % Plot only "visible" sources (with delays <H_length or impulse response peaks >-60dB)
                        h_max_dB60=round(60+20*log10(p_isource(is,ps,1))); %Index of h_max from bottom (1) to top (60) of the colormap
                        if h_max_dB60 <1, h_max_dB60=1; elseif h_max_dB60 >60, h_max_dB60=60; end;
                        mark_colour= cm_LR(h_max_dB60,:);
                        imagesM_h=plot(source(1,is,ps),source(2,is,ps),'+','MarkerEdgeColor',mark_colour); %Plot the visible image source positions (xy,ImageNo,Source)
                        if stem_F==1,stem3(source(1,is,ps),source(2,is,ps),h_max_dB60,'g+' );end; % Plot the visible image source as stem
                    elseif vis_F % Plot "invisible" image sources in black
                        invis_h=plot(source(1,is,ps),source(2,is,ps),'k.'); %Plot positions of invisible image sources (xyz,ImageNo,ParentSource)
                    end;
                    temp_x= source(1,is,ps); temp_y= source(2,is,ps);
                    x_size=[x_size temp_x]; y_size=[y_size temp_y]; % Accumulate size for plotting
                end; %Slice select
            end; % For each image source
        end; % For each parent source
        
    case {'two_mic','mithrir','cipicir'} % Assumes two-sensor inter-sensor axis is parallel with y axis and lies in z=constant plane
        plot(sensor_xyz(1,:),sensor_xyz(2,:),'k'); %Draw line between L & R sensors
        plot(receiver(1),receiver(2),'.k'); % Mark the receiver reference point
        sensorL_h=plot(sensor_xyz(1,1),sensor_xyz(2,1),'o','MarkerEdgeColor','k','MarkerFaceColor','b'); % Position of L sensor     
        sensorR_h=plot(sensor_xyz(1,2),sensor_xyz(2,2),'o','MarkerEdgeColor','k','MarkerFaceColor','r'); % Position of R sensor       
        
        cm_L = L_colormap; % Get the 60 colormap RGB values for left sensor only
        cm_R = R_colormap; % Get the 60 colormap RGB values for right sensor only
        for ps=1:n_sources % For each parent source
            for is=1:n_images % For each image source (with delays <H_length or impulse response peaks >-60dB)
                if ( (source(3,is,ps) > slice_z*(0.99) & source(3,is,ps) < slice_z*(1.01)) ), % Display a slice at the chosen height
                    temp_x=[]; temp_y=[]; %Clear plot size accumulators
                    if p_isource(is,ps,1)>0 % Plot sources "visible" to L sensor
                        h_max_dB60=round(60+20*log10(p_isource(is,ps,1))); %Index from bottom (1) to top (60) of the colormap
                        if h_max_dB60 <1, h_max_dB60=1; elseif h_max_dB60 >60, h_max_dB60=60; end;
                        mark_colour= cm_L(h_max_dB60,:);
                        imagesL_h=plot(source(1,is,ps),source(2,is,ps),'+','MarkerEdgeColor',mark_colour); %Plot the visible image source positions (+) (xy,ImageNo,Source)
                        if stem_F==1,stem3(source(1,is,ps),source(2,is,ps),h_max_dB60,'b+');end; %Plot the visible image source positions (+) (xy,ImageNo,Source)                           
                    end;
                    if p_isource(is,ps,2)>0 % Plot sources "visible" to R sensor
                        h_max_dB60=round(60+20*log10(p_isource(is,ps,2))); %Index from bottom (1) to top (60) of the colormap
                        if h_max_dB60 <1, h_max_dB60=1; elseif h_max_dB60 >60, h_max_dB60=60; end;
                        mark_colour= cm_R(h_max_dB60,:);
                        imagesR_h=plot(source(1,is,ps),source(2,is,ps),'s','MarkerEdgeColor',mark_colour); %Plot the visible image source positions (square) (xy,ImageNo,Source)
                        if stem_F==1,stem3(source(1,is,ps),source(2,is,ps),h_max_dB60,'rs');end; %Plot the visible image source positions (square) (xy,ImageNo,Source)
                    end;
                    if vis_F && (p_isource(is,ps,1)== 0 && p_isource(is,ps,2) == 0) % Plot "invisible" image sources in black
                        invis_h=plot(source(1,is,ps),source(2,is,ps),'k.'); %Plot positions of invisible image sources (xyz,ImageNo,ParentSource)
                    end
                    temp_x= source(1,is,ps); temp_y= source(2,is,ps);
                    x_size=[x_size temp_x]; y_size=[y_size temp_y]; % Accumulate size for plotting
                end; %Slice select
            end; % For each image source
        end; % For each parent source
        
    otherwise
        disp('Unknown receiver set up at freq response plot');
        return
end; % Switch sensor

%Scale the plot to show all image sources and the complete room.
xmin= min([x_size -Lx])-Lx; xmax= max([x_size Lx]+Lx);
ymin= min([y_size -Ly])-Ly; ymax= max([y_size Ly]+Ly);
% Convert max and min to integer number of room lengths
xr_min=floor(xmin/Lx);xr_max=ceil(xmax/Lx);
yr_min=floor(ymin/Ly);yr_max=ceil(ymax/Ly);

% Plot the imaged rooms
for xrooms=xr_min:xr_max
    plot([xrooms xrooms]*Lx,[yr_min yr_max]*Ly,'k:');
end
for yrooms=yr_min:yr_max
    bounds_h=plot([xr_min xr_max]*Lx,[yrooms yrooms]*Ly,'k:');
end;

if stem_F==0, % 2D Plan view
    [rows cols]=size(p_isource);
    if rows*cols > n_sources % If not anechoic room
        rad=H_length*c*T; % Calculate radius corresponding to longest impulse response length
        rectangle('Curvature',[1 1],'Position',[xp-rad,yp-rad,2*rad,2*rad],'EdgeColor','r','LineStyle',':'); %Plot circular limit of audible image sources
    end;
    axis([xr_min*Lx xr_max*Lx yr_min*Ly yr_max*Ly]); % Scale primary axes in metres
    title(['2D Plot of Room, Image Source(s) and Image Room Boundaries at slice height ' num2str(slice_z) ' m']);
    axis equal;
else, % 3D stem plot
    axis([xr_min*Lx xr_max*Lx yr_min*Ly yr_max*Ly 0 60]); % Scale primary axes in metres and height in dB
    title(['3D Stem plot of Room, Image Source(s) and Image Room Boundaries at slice height ' num2str(slice_z) ' m']);
    zlabel('Intensity (NB +60 dB offset)');
    view(3); % Set view for 3D stem plot
end;
xlabel('Length (x) m'); ylabel('Width (y) m');

switch sensor % ---------- Identify the receiver system and display the legend ---------------------
    case {'one_mic'}
        if vis_F, %Show legend that features "invisible" image sources
            legend([outline_h; bounds_h; source_h; sensorM_h; imagesM_h; invis_h],{'Ay1';'Ax1';'Ay2';'Ax2'...
                    ;'Room Bounds';'Source(s)';'Receiver';'Images Mono';'"Inaudible"'},-1);
        else
            legend([outline_h; bounds_h; source_h; sensorM_h; imagesM_h],{'Ay1';'Ax1';'Ay2';'Ax2'...
                    ;'Room Bounds';'Source(s)';'Receiver';'Images Mono'},-1);
        end;
    case {'two_mic','mithrir','cipicir'} % Assumes two-sensor inter-sensor axis is parallel with y axis and lies in z=constant plane
        if vis_F, %Show legend that features "invisible" image sources
            legend([outline_h; bounds_h; source_h; sensorL_h; sensorR_h; imagesL_h; imagesR_h; invis_h],{'Ay1';'Ax1';'Ay2';'Ax2'...
                    ;'Room Bounds';'Source(s)';'Receiver L';'Receiver R';'Images L';'Images R';'"Inaudible"'},-1);
        else
            legend([outline_h; bounds_h; source_h; sensorL_h; sensorR_h; imagesL_h; imagesR_h],{'Ay1';'Ax1';'Ay2';'Ax2'...
                    ;'Room Bounds';'Source(s)';'Receiver L';'Receiver R';'Images L';'Images R'},-1); 
        end;
end;

hold off
drawnow; %Force completion of all previous figure drawing before continuing
%--------------------------- End of roomplot_2D.m -----------------------------------------------
