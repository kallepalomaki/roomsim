function roomsim
%Usage: roomsim;   This Dispatcher presents the welcome notice and the master menu for selecting Roomsim operations
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
%
% Functions called: roomsim_welcome, roomsim_info.m, roomsim_run.m, roomsim_convolve.m, roomsim_disp.m

%*** Fresh Start ***
%(Comment next two lines out to enable breakpoints during debug)
close all; clear all;	%Clear graphics, clear workspace.
warning off; % Suppress warnings

%*** Global Declarations ***
global curr_fig; %Keep track of current figure number
global log_fid; % Identifier of logfile

%*** Constants ***
expiry_day='1 July 2003';
%-----------------------------------------------

roomsim_welcome;% Display Welcome Screen

if datenum(expiry_day) - fix(now) <= 0,
    expiry_message=strvcat(...
        'A new version of this program should now be available'...
        ,'You are advised to replace this BETA version '...
        );
    h = warndlg(expiry_message,'Roomsim');
    beep;
    uiwait(h);% Wait for user to acknowledge
    quit; % Exit from MATLAB
end;

%*** Now we Start ***
log_fid =fopen('Roomsim.log','at'); % Open the logfile for writing text and get an identifier
fprintf(log_fid,['\n\nOpening Roomsim   ' datestr(now) '\n\n']); % Stamp the date and time of opening on the log file

M_banner='Roomsim: Main Menu';
Button1='Exit to MATLAB (NB Clears windows)';
Button2='About Roomsim';
Button3='Roomsim Licence';
Button4='Set-up and Run the Room Simulation';
Button5='Utilities';
Button6='Display Options';
Button7='Vacancy';
Button8='Vacancy';
choice_m0 = 0;
while choice_m0 ~= 1,
    choice_m0=menu(M_banner,Button1,Button2,Button3,Button4,Button5,Button6,Button7,Button8);
    beep;
    switch choice_m0
        case 2, % *** Display programme info ***
            roomsim_info;
            beep;
            
        case 3,	% *** Set-up and run the room simulation ***
            roomsim_licence;
            
        case 4,	% *** Set-up and run the room simulation ***
            roomsim_run;
            beep;
            
        case 5,%*** Utilities ***
            roomsim_utils;
            beep;
            
        case 6,%*** Display Options ***
            roomsim_disp;
            beep;
            
        case 7,%*** Expansion space ***
            h=msgbox('This item has been left vacant for future expansion','Roomsim');
            beep;
            uiwait(h);% Wait for user to acknowledge
            
        case 8, %*** Expansion space ***
            h=msgbox('This item has been left vacant for future expansion','Roomsim');
            beep;
            uiwait(h);% Wait for user to acknowledge
    end;%of menu0 choices
    
end; % of menu0 while loop

%*** Now clean up and shut down ***
fprintf(log_fid,['\n\nClosing Roomsim   ' datestr(now) '\n\n']); % Stamp the date and time of closing on the log file
fclose(log_fid); % Close the log file
colormap(JET); % Reset colour lookup table
close all; % Shut down all the windows
clear all; % Clear the workspace			

%-------------------- End of roomsim.m ------------------------------
