function roomsim_utils
%Usage: roomsim_utils;   This Dispatcher allows the user to select utility functions within Roomsim
%------------------------------------------------------------------
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
%-------------------------------------------------------------------------
% Functions called: roomsim_high_pass.m, roomsim_low_pass.m, roomsim_convolve.m, calc_RT60.m

%*** Global Declarations ***
global curr_fig; %Keep track of current figure number

M_banner='Utilities Menu';
Button1='Exit from Utilities';
Button2='Remove DC & LF ripple (High Pass Filter Cut-off = 0.004*Nyquist Frequency)';
Button3='Anti-alias (Low Pass Filter Cut-off = 0.9*Nyquist Frequency)';
Button4='Convolve Impulse Response & Audio Data';
Button5='Convert Roomsim audio files *.mat to *.wav or vice versa';
Button6='Vacancy';
Button7='Vacancy';
Button8='Vacancy';
choice_m2 = 0;
while choice_m2 ~= 1,
    choice_m2=menu(M_banner,Button1,Button2,Button3,Button4,Button5,Button6,Button7,Button8);
    beep;
    
    fig_loc=[50 50 450 350]; % Location for following plots
    switch choice_m2
        case 2, % *** High-Pass filter used for LF ripple & DC removal ***
            high_pass_filter(fig_loc); % High-Pass 4th order Butterworth filter Cut-off 0.0023*Nyquist frequency
            beep;
            
        case 3,	% ***  Low-Pass used for bandlimiting/anti-aliasing ***
            low_pass_filter(fig_loc); % Low-Pass 4th order Butterworth filter Cut-off 0.91*Nyquist frequency
            beep;
            
        case 4,	% *** Convolve two files ***
            roomsim_convolve(fig_loc);
            beep;
            
        case 5,%*** Convert from *.mat to *.wav or vice versa ***
            roomsim_convert;
            beep;
            
        case 6,%*** Expansion space ***
            h=msgbox('This item has been left vacant for future expansion','Roomsim Utilities');
            beep;
            uiwait(h);% Wait for user to acknowledge
            
        case 7,%*** Expansion space ***
            h=msgbox('This item has been left vacant for future expansion','Roomsim Utilities');
            beep;
            uiwait(h);% Wait for user to acknowledge
            
        case 8, %*** Expansion space ***
            h=msgbox('This item has been left vacant for future expansion','Roomsim');
            beep;
            uiwait(h);% Wait for user to acknowledge
    end;%of menu2 choices
    
end; % of menu2 while loop
