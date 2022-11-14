function [Fs,c,humidity,order,H_length,H_filename,air_F,smooth_F,Fc_HP,plot_F2,plot_F3,dist_F,alpha_F...
        ,room_size,receiver,sensor,sensor_space,sensor_dir,S_No,source_polar,F_abs,A]=roomsim_setup
%Usage: [Fs,c,humidity,order,H_length,H_filename,air_F,smooth_F,Fc_HP,plot_F2,plot_F3,dist_F,alpha_F...
%         ,room_size,receiver,sensor,sensor_space,sensor_dir,S_No,source_polar,F_abs,A]=roomsim_setup;
% Setup code for Roomsim. This function prompts for either:
% 1) A previously saved *.mat file with data entered through the Manual Input dialogue
% 2) Input from an Excel spreadsheet file (text.xls can be edited by the user)
% 3) Manual input via dialogue prompts 
% Choices 2 and 3 above will be given the option to save resulting setup data as a MAT file for future use.
%------------------------------------------------------------------------------------- 
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
%----------------------------------------------------------------------------------------
%
% The required parameters for this function are:
% Fs, the sampling frequency (Hz).
% humidity, of the air in the room. A formula (valid for 20<= h <= 70) is used to compute the absorption m of the air volume (%).
% temperature, of the air in the room. This is used to compute the speed of sound in the room (m/s)
% order, a column vector containing the order of reflections to be calculated for each dimension x,y,z.
% H_length, maximum required length of room impulse response(s).
% H_filename, a string identifier naming the file to which the final impulse responses H will be saved.
% air_F, logical Flag value 0 or 1, false = no absorption due to air, true = air absorption is present.
% smooth_F, logical Flag value 0 or 1, false = no smoothing filter applied to weight intersample impulses, true = smoothing filter used (A5)
% Fc_HP, scalar, the cut-off frequency of the high-pass filter applied to attenuate DC, 0 = no High-pass filter (A6 used 100 Hz)
% plot_F2, logical Flag value 0 or 1, Set this to plot the geometry of receiver, source(s), image sources & image rooms,
%   false= no plot, true= 2D-plan showing image rooms on constant z plane.
% plot_F3, logical Flag value 0 or 1, Set for a rotatable 3D-plot using the MATLAB controls.
% dist_F, logical Flag value 0 or 1, Set to enable 1/R distance attenuation effect.
% alpha_F, logical Flag value 0 or 1, Set to enable display of room surface opacity proportional to reflectivity
% room_size(Lx,Ly,Lz), a column vector containing the length x width y and height z of the room in metres.
% receiver(xp,yp,zp), a column vector containing the RH Cartesian cordinates locating the receiver head in metres.
%   (NB no facility to set an angle for the Head sagittal plane ie it is aligned with the x axis).
% sensor, text string value identifying receiver system, 
%   one of: One sensor ('one_mic'), Two sensors (two_mic'), MIT Kemar ('mithrir'), CIPIC Head ('cipicir').
% Microphone directionality (deg) (Omni is -180,180,-90,90)
%   min_azim_sensor, a scalar, minimum azimuth (deg) viewable by the sensor
%   max_azim_sensor, a scalar, maximum azimuth (deg) viewable by the sensor
%   min_elev_sensor, a scalar, minimum elevation (deg) viewable by the sensor
%   max_elev_sensor, a scalar, maximum elevation (deg) viewable by the sensor
% sensor_space, a scalar, sensor separation for the two_sensor case (if MIT or CIPIC it is implicit in the HRIR data)
% S_No, text value, CIPIC subject number format '&&&' (e.g. '021' is Kemar with small pinnae)
% R_s, radial distance (m) between receiver and source(s). NB a check is performed to ensure that all sources are inside the room.
% alpha, a scalar (or vector) of azimuth(s), -180< alpha < 180 (deg), used to set the source location(s) relative to receiver.
%   NB 0 deg has a source located in the xz plane , +ve deg rotates Anti-CW on the xy plane viewed in plan.
% beta, a scalar (or vector) of elevations(s), -90< beta < 90 (deg), used to set the source location(s) relative to receiver.
%   NB 0 deg has a source located in the xy plane, +ve deg rotates upwards.
% A, a matrix of Sabine energy surface absorption coefficients for the "shoebox" enclosure containing (in columns)
%   Ax1,Ax2,Ay1,Ay2,Az1,Az2, all being column vectors of frequency dependent absorption coefficients at the six standard measurement frequencies.
%   x1, y1 and z1 refer to the surfaces lying in the x=0, y=0 and z=0 planes respectively, while
%   x2, y2 and z2 refer to the surfaces lying in the x=Lx, y=Ly and z=Lz planes respectively.
%   Thus x1 & x2, y1 & y2, z1 & z2 are opposing surfaces eg. z1 refers to the floor and z2 the ceiling.
%-----------------------------------------------------------------------------------
% Functions called: excel_setup.m, 

%*** Globals ***
global log_fid; % Identifier of logfile
global hrir_l hrir_r; % Globals for CIPIC impulse responses left and right
global H3D; % Global for MIT impulse response data

% Declarations
S_No=[]; % Create array for CIPIC subject number
min_elev_sensor=[]; max_elev_sensor=[]; min_azim_sensor=[]; max_azim_sensor=[]; % Create arrays for sensor directionality.
sensor_space=[]; % Create array for sensor separation
error_title='Roomsim_setup: Error';

%Display a menu
M_banner='Roomsim_setup: Set up the simulation parameter values';
Button1='Load previously saved MAT file';
Button2='Read from Excel spreadsheet';
Button3='Manual input following prompts';

Fs = 0;
while Fs==0, % Loop until data has been entered in some manner
    M_Inputcase = menu(M_banner,Button1,Button2,Button3);
    beep;
    switch M_Inputcase
        
        case 1 %Input from MAT file
            [name path] = uigetfile('SAVE_roomsim.mat', 'Roomsim_setup: Select a MAT File containing set-up data',0,0); %Display the dialogue box
            if ~any(name),
                %**Alternate return for cancel operation
            else
                load([path name]); %Put the previously stored values into the workspace
                return; ; % Nothing more to do so exit to calling program
            end
            
        case 2 %Input from Excel spreadsheet file setup.xls
            
            [name path] = uigetfile('setup.xls', 'Roomsim_setup: Select an Excel File containing set-up data',0,0); %Display the dialogue box
            if ~any(name), 
                %**Alternate return for cancel operation
            else
                [Fs,humidity,temperature,order,H_length,H_filename,air_F,smooth_F,Fc_HP,plot_F2,plot_F3,dist_F,alpha_F...
                        ,c,room_size,receiver,sensor,sensor_space,sensor_dir,S_No,source_polar,F_abs,A]=excel_setup([path name]);
            end
            
            
        case 3 %Manual input
            %----------Start of user input section-----------------------------
            answer={};
            while isempty(answer),
                banner = 'Roomsim_setup: Enter Simulation control parameters';
                prompt = {'Sampling frequency Fs > 8000 (Hz) :'...
                        ,'Humidity of air 20<= h <= 70 (%):'...
                        ,'Temperature of air (Celcius):'...
                        ,'Limit to Order of reflections (-1 Program decides):'...
                        ,'Limit to Impulse response length (-1 Program decides):'...
                        ,'Filename for Impulse response:'...
                        ,'Air flag, 1 = Air absorption present (0 = not present):'...
                        ,'Smoother flag, 1 = Smoother present (0 = not present):'...
                        ,'High-Pass filter cut-off (Hz), scalar value eg 50~100 (0 = filter not present):'...
                        ,'2D Plotting flag 1 = Display 2D Plot (0 = No Plot):'...
                        ,'3D Plotting flag 1 = Display 3D Plot (0 = No Plot):'...
                        ,'Distance flag, 1 = Distance Attenuation present (0 = not present):'...
                        ,'Transparency flag, 1 = not opaque (0 = reflectivity sets opacity):'...
                    };
                lines = 1;
                def = {'44100','50','20','-1','-1','ROOM_IMPULSE','0','0','0','0','0','1','0'}; %Default values
                answer = inputdlg(prompt,banner,lines,def);
                beep; 
            end; % Of while loop to disable inappropriate CANCEL button operation
            Fs=str2num(answer{1});
            humidity=str2num(answer{2});
            temperature=str2num(answer{3});
            order=str2num(answer{4});
            H_length=str2num(answer{5});
            H_filename=answer{6};
            air_F=logical(str2num(answer{7}));
            smooth_F=logical(str2num(answer{8}));
            Fc_HP=logical(str2num(answer{9}));
            plot_F2=logical(str2num(answer{10}));
            plot_F3=logical(str2num(answer{11}));
            dist_F=logical(str2num(answer{12}));
            alpha_F=logical(str2num(answer{13}));
           
            c=round(331*sqrt(1+0.0036*temperature)); % Calculate speed of sound (m/s) for use in sensor separation check
            
            % ---------------------------- Room, receiver and sources details ----------------------
            % Room dimensions
            Lz=1.25; % Basis length (m) for default values 
            def_Lx=num2str(5*Lz);
            def_Ly=num2str(3*Lz);
            def_Lz=num2str(2*Lz);
            def_vol=num2str(round(5*3*2*Lz^3)); % Compute volume of default room
            
            answer={};
            while isempty(answer),
                banner = ['Room size (Default is ' def_vol 'm^3 room, Volkmann 5:3:2 ratio)'];
                prompt = {'Enter Length (Depth) (Lx) of room in meters (m) :'...
                        ,'Enter Width (Ly) of room in meters (m) :'...
                        ,'Enter Height (Lz) of room in meters (m) :'};
                lines = 1;
                def = {def_Lx,def_Ly,def_Lz}; %Default values (Volkmann 2:3:5 ratio 'living room/3 Person office dimensions')
                answer = inputdlg(prompt,banner,lines,def,'on');
                beep;
            end; % Of while loop to disable inappropriate CANCEL button operation
            Lx=str2num(answer{1});
            Ly=str2num(answer{2});
            Lz=str2num(answer{3});
            
            % Receiver (sensor) details
            answer={};
            while isempty(answer),
                banner = 'Receiver reference position [xp,yp,zp](m)';
                prompt = {'Enter Receiver x co-ordinate (m) (Default is room_length/4) :'...
                        ,'Enter Receiver y co-ordinate (m) (Default is room_width/2) :'...
                        ,'Enter Receiver z co-ordinate (m) (Default is typical seated ear height) :'};
                lines = 1;
                def = {num2str(Lx/4),num2str(Ly/2),'1.1'}; %Default values
                answer = inputdlg(prompt,banner,lines,def);
                beep;
            end; % Of while loop to disable inappropriate CANCEL button operation
            xp=str2num(answer{1});
            yp=str2num(answer{2});
            zp=str2num(answer{3});
            
            ButtonName=0;
            while ButtonName==0, % Loop until valid choice has been made (disable close window button)
                ButtonName=menu('Specify Receiver System', ...
                    'One sensor','Two sensors','MIT Kemar','CIPIC Head');
                beep;
            end; % Of while loop to disable inappropriate close window button operation
            switch ButtonName,
                case 1, 
                    sensor='one_mic';
                    out_of_range = true;
                    while out_of_range, %Loop if sensor directionality is not sensible
                        answer={};
                        while isempty(answer),
                            banner = 'Sensor Directional Characteristics (Default is Omni-directional)';
                            prompt = {'Enter minimum azimuth :'...
                                    ,'Enter maximum azimuth :'...
                                    ,'Enter minimum elevation :'...
                                    ,'Enter maximum elevation :'};
                            lines = 1;
                            def = {'-180','180','-90','90'}; %Default values for Omnidirectional
                            answer = inputdlg(prompt,banner,lines,def);
                            beep;
                        end; % Of while loop to disable inappropriate CANCEL button operation
                        min_azim_sensor=str2num(answer{1});
                        max_azim_sensor=str2num(answer{2});
                        min_elev_sensor=str2num(answer{3});
                        max_elev_sensor=str2num(answer{4});
                        % Check values within range
                        if min_elev_sensor < -90 || max_elev_sensor > 90
                            message='-90< Sensor elevation <90 deg';
                        elseif min_elev_sensor >= max_elev_sensor
                            message='Minimum elevation must be less than Maximum';
                        elseif min_azim_sensor < -180 || max_azim_sensor > 180
                            message='-180< Sensor azimuth <180 deg';
                        elseif min_azim_sensor >= max_azim_sensor
                            message='Minimum azimuth must be less than Maximum';
                        else
                            out_of_range = false; % Input accepted
                        end;
                        if out_of_range  
                            banner='Re-enter value';
                            h=msgbox(message,banner,'warn');  %Warn & beep.
                            beep;
                            uiwait(h);% Wait for user to acknowledge
                            answer={};
                        end;
                    end;
                    
                case 2,
                    sensor='two_mic';
                    delay_m=0;
                    while delay_m <1, %Loop if sensors are separated by less than one sample delay
                        answer={};
                        while isempty(answer),
                            banner = 'Two sensor array';
                            prompt = {'Enter Sensor separation (m) (Default is CIPIC average head width) :'};
                            lines = 1;
                            def = {'0.145'}; %Default value is CIPIC average head width rounded to nearest mm.
                            answer = inputdlg(prompt,banner,lines,def);
                            beep;
                        end; % Of while loop to disable inappropriate CANCEL button operation
                        sensor_space=str2num(answer{1});
                        delay_m=sensor_space*Fs/c; % delay in samples = (Samples per sec)*(Distance)/speed of sound
                        if delay_m < 1, % Check two sensor separation is greater than one sample distance
                            banner='Re-enter value';
                            message='Two sensor separation < 1 sample distance';
                            h=msgbox(message,banner,'warn');  %Warn & beep.
                            beep;
                            uiwait(h);% Wait for user to acknowledge
                        end;
                    end;
                    
                    out_of_range = true;
                    while out_of_range, %Loop if sensor directionality is not sensible
                        answer={};
                        while isempty(answer),
                            banner = 'Sensor Directional Characteristics (Default is Omni-directional)';
                            prompt = {'Enter minimum azimuth [L sensor, R sensor] :'...
                                    ,'Enter maximum azimuth [L sensor, R sensor] :'...
                                    ,'Enter minimum elevation [L sensor, R sensor] :'...
                                    ,'Enter maximum elevation [L sensor, R sensor] :'};
                            lines = 1;
                            def = {'[-180 -180]','[180 180]','[-90 -90]','[90 90]'}; %Default values for Omnidirectional
                            answer = inputdlg(prompt,banner,lines,def);
                            beep;
                        end; % Of while loop to disable inappropriate CANCEL button operation
                        min_azim_sensor=str2num(answer{1});
                        max_azim_sensor=str2num(answer{2});
                        min_elev_sensor=str2num(answer{3});
                        max_elev_sensor=str2num(answer{4});
                        % Check values within range
                        if min(min_elev_sensor) < -90 || max(max_elev_sensor) > 90
                            message='-90< Sensor elevation <90 deg';
                        elseif any(min_elev_sensor >= max_elev_sensor)
                            message='Minimum elevation must be less than Maximum';
                        elseif min(min_azim_sensor) < -180 || max(max_azim_sensor) > 180
                            message='-180< Sensor azimuth <180 deg';
                        elseif any(min_azim_sensor >= max_azim_sensor)
                            message='Minimum azimuth must be less than Maximum';
                        else
                            out_of_range = false; % Input accepted
                        end;
                        if out_of_range, %Warn & beep and loop back for re-entry of data.
                            banner='Re-enter value';
                            h=msgbox(message,banner,'warn');  
                            beep;
                            uiwait(h);% Wait for user to acknowledge
                            answer={};
                        end;
                    end;
                case 3,
                    sensor='mithrir';
                    %---------- Extract hrir pairs from MIT data base and load into workspace cell array H3D --------
                    MIT_root = 'MIT_HRTF'; %Root name of MIT directory (folder)
                    if (exist(MIT_root,'dir')==7) % Check path name to the MIT Kemar directory (folder)
                        H3D={}; % Declare array for loading MIT data
                        
                        dir_ch = '\'; 
                        subdir='Kemar\compact';
                        filename='hrir_final';
                        ext = '.mat';
                        MIT_file = [MIT_root dir_ch subdir dir_ch filename ext]; %Form pathname to Kemar data file
                        
                        load(MIT_file); % Load the MIT Kemar data file into the workspace (creates H3D)
                    else
                        h=errordlg([MIT_file '  NOT FOUND, exiting from roomsim'],error_title);
                        beep;
                        uiwait(h);
                        return;
                    end
                    
                case 4,
                    sensor='cipicir';
                    CIPIC_root = 'CIPIC_HRTF'; %Root name of CIPIC directory (folder)
                    if (exist(CIPIC_root,'dir')==7)
                        S_No=-1;
                        while S_No < 0,
                            ITD=[]; OnL=[]; OnR=[]; hrir_l=[]; hrir_r=[]; % Declare arrays for loading CIPIC data
                            ok=0;
                            while ok==0,
                                % Setup list box for picking CIPIC subject number           
                                subject_list={'003','008','009','010','011','012','015','017','018','019','020'...
                                        ,'021','027','028','033','040','044','048','050','051','058','059'...
                                        ,'060','061','065','119','124','126','127','131','133','134','135','137'...
                                        ,'147','148','152','153','154','155','156','158','162','163','165'};
                                
                                [selection,ok] = listdlg('PromptString','Select a CIPIC subject number (021 is Kemar small pinnae)'...
                                    ,'SelectionMode','single','InitialValue',12,'ListSize',[300 300],'ListString',subject_list);
                                beep;
                            end; % Of while loop to disable inappropriate CANCEL button operation
                            S_No = subject_list{selection};% Return CIPIC subject number (as character)
                            
                            dir_ch = '\'; 
                            subdir=['standard_hrir_database\subject_' S_No];
                            filename='hrir_final';
                            ext = '.mat';
                            CIPIC_file = [CIPIC_root dir_ch subdir dir_ch filename ext]; %Form pathname to subject file
                            if (exist(CIPIC_file,'file')==2), % Check selected subject file is installed
                                load(CIPIC_file); % Load the subject file into the workspace (creates hrir_l and hrir_r)
                            else
                                h=errordlg([CIPIC_file '  NOT FOUND, check Subject No ' S_No 'is installed.'],error_title);
                                beep;
                                uiwait(h);
                                S_No=-1;
                            end;
                        end; % Of while loop to catch missing subject number
                    else
                        h=errordlg([CIPIC_root '  NOT FOUND, exiting from roomsim'],error_title);
                        beep;
                        uiwait(h);
                        return;
                    end
                    
            end % receiver switch
            
            % Source location(s)
            answer={};
            while isempty(answer),
                banner = 'Source position(s) relative to Receiver (Polar Co-ordinates)';
                prompt = {'Enter Source radial distance(s) [R1 R2 R3 ..] (m) :'...
                        ,'Enter Source azimuth(s) [az1 az2 az3 ..] (-180<deg<180) :'...
                        ,'Enter Source elevation(s) [el1 el2 el3 ..] (-90<deg<90) :'};
                lines = 1;
                def = {'[1]','[0]','[0]'}; %Default values
                answer = inputdlg(prompt,banner,lines,def);
                beep;
            end; % Of while loop to disable inappropriate CANCEL button operation
            R_s=str2num(answer{1});
            alpha=str2num(answer{2});
            beta=str2num(answer{3});
            
            %--------------- Get the absorption coefficients for each surface (and check they exist) ---------
            Abs=[]; Ax1=[]; Ax2=[]; Ay1=[]; Ay2=[]; Az1=[]; Az2=[]; %Initialise absorption coefficient vectors
            surf_error=strvcat('Surface type not recognised','Correct your choice (see surfaces.m)');
            button=[];
            while isempty(button) % Loop until valid choice has been made (disable close window button)
                button = questdlg('All surfaces have identical absorption?','Roomsim_setup:','Yes','No','Yes');
                beep;
            end; % Of while loop to disable inappropriate close window button operation
            switch button
                case 'Yes' % All surfaces identical
                    while isempty(Abs),
                        [Abs F_abs]=get_abs_table('Abs');
                        if isempty(Abs),
                            h=errordlg(surf_error,error_title);
                            beep;
                            uiwait(h);
                        end;
                    end;
                    Ax1=Abs; Ax2=Abs; Ay1=Abs; Ay2=Abs; Az1=Abs; Az2=Abs; 
                case 'No' % Non-identical absorption, get each one individually
                    while isempty(Ax1),
                        [Ax1 F_abs]=get_abs_table('Ax1');
                        if isempty(Ax1),
                            h=errordlg(surf_error,error_title);
                            beep;
                            uiwait(h);
                        end;
                    end;
                    while isempty(Ax2),
                        [Ax2 F_abs]=get_abs_table('Ax2');
                        if isempty(Ax2),
                            h=errordlg(surf_error,error_title);
                            beep;
                            uiwait(h);
                        end;
                    end;
                    while isempty(Ay1),
                        [Ay1 F_abs]=get_abs_table('Ay1');
                        if isempty(Ay1),
                            h=errordlg(surf_error,error_title);
                            beep;
                            uiwait(h);
                        end;
                    end;
                    while isempty(Ay2),
                        [Ay2 F_abs]=get_abs_table('Ay2');
                        if isempty(Ay2),
                            h=errordlg(surf_error,error_title);
                            beep;
                            uiwait(h);
                        end;
                    end;
                    while isempty(Az1),
                        [Az1 F_abs]=get_abs_table('Az1');
                        if isempty(Az1),
                            h=errordlg(surf_error,error_title);
                            beep;
                            uiwait(h);
                        end;
                    end;
                    while isempty(Az2),
                        [Az2 F_abs]=get_abs_table('Az2');
                        if isempty(Az2),
                            h=errordlg(surf_error,error_title);
                            beep;
                            uiwait(h);
                        end;
                    end;
            end;
            
            %---------------- End of surface absorption manual set up ---------------
            % Pack various parameters for the room simulation
            room_size=[Lx;Ly;Lz]; % Pack up room dimensions into column vector.
            receiver=[xp;yp;zp]; % Pack up receiver (listener's head) coordinates into column vector.
            sensor_dir=[min_azim_sensor; max_azim_sensor; min_elev_sensor; max_elev_sensor]; % Pack up sensor directionality one column per sensor
            source_polar=[R_s;alpha;beta]; % Pack up source(s) coordinates into column vector.
            A=[Ax1 Ax2 Ay1 Ay2 Az1 Az2]; % Pack up column vectors of absorption coefficients in array A
    end; % of M_inputcase switch
end; % of menu while loop

%----------------- Save set-up data to a MATLAB loadable file-----------------
beep;
file_spec='SAVE_Roomsim.mat';
title='Roomsim_setup: Save data as a MAT file';
dot_ext=[];
while isempty(dot_ext),
    [name path] = uiputfile(file_spec, title); %Display the dialogue box
    if ~any(name), % File select was cancelled, exit the while loop.
        fprintf(log_fid,'\n\n Roomsim_setup: Setup parameters have not been saved'); % Print to log file (avoids wait for user response)
        break; % Alternate return for cancel operation. No data saved to file. 
    end;
    filename = lower([path name]); % Force path+filename to lowercase
    dot_ext = findstr(filename,'.'); % Find the extension marker
    if isempty(dot_ext) | length(filename)~=dot_ext+3 | ~strcmp(filename(dot_ext+1:dot_ext+3),'mat'),
        h=errordlg('You must specify a 3 letter *.mat extension','roomsim_run error');
        beep;
        uiwait(h);
        dot_ext=[]; % Force resubmit of filename
    else,
        switch filename(dot_ext+1:dot_ext+3),
            case 'mat'
                % Using the functional form of SAVE (NB Non-EVAL code for compilation)
                save([path name]); % Save all variables in workspace to a ".mat" file
                fprintf(log_fid,'\n\n Roomsim_setup: Setup parameters have been saved to %s', [path name]); % Print to log file (avoids wait for user response)
            otherwise,
                h=errordlg('Data file extension not recognised, exiting','savef error');
                beep;
                uiwait(h);
                return;
        end;
    end;
end; % Of while loop to trap missing .mat extension
%-------------------------------- End of roomsim_setup.m ----------------------------------------
