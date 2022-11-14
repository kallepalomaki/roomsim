function roomsim_disp
% Usage: roomsim_disp.m; This function allows user to plot saved roomsim data without re-running roomsim_core.
%--------------------------------------------------------------------------- 
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
%-------------------------------------------------------------------------------
%Functions called: roomplot_imp.m, roomplot_magf.m, roomplot_2D.m, roomplot_3D.m

%*** Global Declarations ***
global curr_fig; %Keep track of current figure number

M_banner='Display Menu: Choose a plot format';
Button1='Exit from Display Options';
Button2='Plot of the room impulse response';
Button3='Plot of the room frequency response';
Button4='2D Plot of the room and the image rooms';
Button5='3D Plot of the room and the image sources';

M_Inputcase =0;
while M_Inputcase ~= 1 % Get the data file and display the plotting menu
    M_Inputcase = menu(M_banner,Button1,Button2,Button3,Button4,Button5);
    beep;        
    if M_Inputcase ~=1,
        [name path] = uigetfile('PLOT_roomsim.mat', 'Select previously saved data for plotting'); %Display the dialogue box
        if ~any(name), 
            return; %**Alternate return for cancel operation
        else
            load([path name]); %Put the previously stored values into the workspace
        end;
    end;
    
    fig_loc=[50 50 450 350]; % Location for following plots
    switch M_Inputcase
        case 2 %-----------------------Plot the impulse responses to each sensor---------------
            roomplot_imp(Fs,H,sensor,fig_loc);
            
        case 3 %-----------------------Plot the frequency responses to each sensor---------------
            roomplot_magf(Fs,H,HRTF,sensor,fig_loc);
            
        case 4 %----------------------- 2D Plot of the room and the image rooms---------------
            % Display the receiver, source, image sources and image rooms as a 2D plan.
            roomplot_2D(c,Fs,room_size,source_xyz,receiver,sensor,sensor_xyz,sensor_dir,H_length,source,p_isource,fig_loc);
            
        case 5 %----------------------- 3D Plot of the room and the image sources---------------
            roomplot_3D(room_size,source_xyz,receiver,sensor,sensor_xyz,sensor_dir,source,p_isource,fig_loc);
            
    end; % of menu choices
    
end; % of menu while loop
%----------------------------------------- End of roomsim_disp.m --------------------------------