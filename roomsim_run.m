function roomsim_run
% Usage: roomsim_run;    Calls the set-up function, does some checks, calls the core image calculation
% function, calls various ploting functions, and allows saving the set-up as a MAT file for later re-use.
%-------------------------------------------------------------------------------- 
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
% Functions called: roomsim_setup.m, roomplot.m, Norris_Eyring.m, roomsim_core.m, roomplot_imp.m,
%   roomplot_magf.m, roomplot_2D.m, roomplot_3D.m .

%*** Global Declarations ***
global curr_fig; %Keep track of current figure number
global log_fid; % Identifier of logfile
global hrir_l hrir_r; % Globals for CIPIC impulse responses left and right
global H3D; % Global for MIT impulse response data

%------------------- Constants -----------------------------
degtorad=pi/180; % Conversion factor degrees to radians
radtodeg=180/pi; % Conversion factor radians to degrees
head_width=0.145; % %Value is CIPIC average head width rounded to nearest mm.

% Common message texts for roomsim_run 
terminate_msg='Roomsim_run: Program terminating'; % Termination title for error dialogue box
error_title='Roomsim_setup: Error';

%------------------------------ Start of user input section -----------------------------
% Call the function to get the user setup values
[Fs,c,humidity,order,H_length,H_filename,air_F,smooth_F,Fc_HP,plot_F2,plot_F3,dist_F,alpha_F...
        ,room_size,receiver,sensor,sensor_space,sensor_dir,S_No,source_polar,F_abs,A]=roomsim_setup;

%------------------------------- End of user input section ---------------------

% Simulation constants dependent on setup
T=1/Fs; %Sampling period
nyquist=Fs/2;
Fs_c = Fs/c; % Samples per metre

% Other initial declarations
fig_loc=[50 50 450 350]; % Location for following plots

%------------- Display the Frequency variation of the surface absorption coefficients --------------
roomplot_absorption(F_abs,A,fig_loc);

%---------------------------------------- Start of sanity checks ------------------------------------------------

%------------------- Check sampling frequency -----------------------
% Check and prevent possible aliasing of HRIR data
switch sensor % ---------- Identify receiver system---------------------
    case {'mithrir','cipicir'}
        if (Fs ~= 44100) % Warn of inappropriate sampling frequency
            msg_title='Roomsim_run: Forcing Fs= 44.1kHz';
            message=['MIT and CIPIC HRIR data requires Fs = 44100 Hz, Your Fs = ' num2str(Fs) ' Hz'];
            h=msgbox(message,msg_title,'warn');  %Warn & beep.
            beep;
            uiwait(h);% Wait for user to acknowledge
            Fs=44100;
        end;
end;

% Check and warn of possible aliasing of absorption frequencies

if max(F_abs)>= nyquist % Warn Fs too low for absorption data
    answer={};
    while isempty(answer),msg_title = 'Roomsim_run: Low sampling frequency';
        prompt = {['Set a new sampling frequency > 2*' num2str(max(F_abs))]};
        lines = 1;
        def = {num2str(2*max(F_abs))}; %Default value
        answer = inputdlg(prompt,msg_title,lines,def);
        beep;
    end; % Of while loop to disable inappropriate CANCEL button operation
    Fs=str2num(answer{1});
end

%--- Convert source radial distance (R_s), azimuth (alpha) and elevation (beta) to x,y,z coordinates of each parent source -------
xp=receiver(1);
yp=receiver(2);
zp=receiver(3);
hypxy=source_polar(1,:).*cos(degtorad*source_polar(3,:)); %projection of R_s on xy plane
x=xp+hypxy.*cos(degtorad.*source_polar(2,:));% sound source x position.
y=yp+hypxy.*sin(degtorad.*source_polar(2,:));% sound source y position.
z=zp+source_polar(1,:).*sin(degtorad*source_polar(3,:));% sound source z position.
source_xyz=[x;y;z]; % Pack up source(s) coordinates into array one column vector per source.

%-------- Check for primary source/receiver coincident ---------------------------
n_sources=size(source_xyz,2); % Number of sources = number of columns in source_xyz

for ps=1:n_sources % For each parent source
    switch sensor, % Select receiver system
        case {'one_mic','mithrir','cipicir'}
            % Calculate distance from Receiver reference point (Centre of Head, single sensor or mid-point of Two sensors ) to each parent source.
            dist_s=norm((source_xyz(:,ps)-receiver(:)),2); %Pythagoras in 3D
            % calculate delay in samples between receiver and source.
            delay_s=Fs_c*dist_s; % delay in samples = (Samples per sec)*(Distance)/speed of sound
        case 'two_mic', % This code segment allows for sources to exist between a widely spaced sensor pair
            %Microphone pair centered on (xp,yp,zp), lying in xy plane and aligned with y axis, separated by a spacing of sensor_space metres
            sensor_xyz(:,1) = receiver + [0;sensor_space/2;0]; % Add half sensor spacing to yp to give coords of L sensor
            sensor_xyz(:,2) = receiver + [0;-sensor_space/2;0]; % Subtract half sensor spacing from yp to give coords of R sensor
            dist_m(1) = norm(source_xyz(:,ps)-sensor_xyz(:,1),2); %Distance (m) between source and receiver left sensor
            dist_m(2) = norm(source_xyz(:,ps)-sensor_xyz(:,2),2); %Distance (m) between source and receiver right sensor
            delay_s = Fs_c*dist_m; % L & R delays in samples = (Samples per metre)*Distance
    end; % switch sensor
    
    % Warn of sources within one sample period of receiver/sensor location
    if ~all(delay_s >= 1) % Some elements of the test are zero so display a menu & beep.
        beep;      
        M_Title=['Roomsim_run: WARNING Source ' num2str(ps) ' and sensor coincident.'];
        B1_text='Continue with present value';
        B2_text='Return to Roomsim Main Menu';
        beep;
        M_IRLcase = menu(M_Title,B1_text,B2_text);
        switch M_IRLcase
            case 1
%                 Continue
            case 2
                return; % Return to calling program (main menu)
        end;
    end;
end;

%---------- Check for sources outside of room and find max and min axis values for plotting ------------
source_outside=false; %Initialise flag
Lx=room_size(1);
Ly=room_size(2);
Lz=room_size(3);

source_error='Source outside of room, '; % Message warning of impossible source location
if min(x)<0
    xmin=min(x);
    h=errordlg([source_error ' x = ' num2str(xmin) '.'],terminate_msg);% Warn of impossible source location
    beep;
    uiwait(h);
    source_outside=true;
else
    xmin=-0.1;
end
if min(y)<0
    ymin=min(y);
    h=errordlg([source_error ' y = ' num2str(ymin) '.'],terminate_msg)% Warn of impossible source location
    beep;
    uiwait(h);
    source_outside=true;
else
    ymin=-0.1;
end
if min(z)<0
    zmin=min(z);
    h=errordlg([source_error ' z = ' num2str(zmin) '.'],terminate_msg);% Warn of impossible source location
    beep;
    uiwait(h);
    source_outside=true;
else
    zmin=-0.1;
end
if max(x)>=Lx
    xmax=max(x);
    h=errordlg([source_error ' x = ' num2str(xmax) '.'],terminate_msg);% Warn of impossible source location
    beep;
    uiwait(h);
    source_outside=true;
else
    xmax=Lx+0.1;
end
if max(y)>=Ly
    ymax=max(y);
    h=errordlg([source_error ' y = ' num2str(ymax) '.'],terminate_msg);% Warn of impossible source location
    beep;
    uiwait(h);
    source_outside=true;
else
    ymax=Ly+0.1;
end
if max(z)>=Lz
    zmax=max(z);
    h=errordlg([source_error ' z = ' num2str(zmax) '.'],terminate_msg);% Warn of impossible source location
    beep;
    uiwait(h);
    source_outside=true;
else
    zmax=Lz+0.1;
end
%------------------- End of source location checks -----------------------

%-------------------- Check for sensors outside of room --------------------
sensor_outside=false; %Initialise flag
switch sensor % ---------- Identify receiver system---------------------
    case {'one_mic'}
        sensor_xyz = receiver; % Single sensor at receiver reference point
        if (min(min(sensor_xyz))<0)||(max(sensor_xyz(1))>=Lx)||(max(sensor_xyz(2))>=Ly)||(max(sensor_xyz(3))>=Lz)
            sensor_error=strvcat('Sensor outside of room.');
            h=errordlg(sensor_error,terminate_msg);% Warn of impossible sensor location
            beep;
            uiwait(h); % Wait for user response
            sensor_outside=true;
        end;
        
    case {'two_mic'} %(Assumes two-sensor inter-sensor axis is parallel with y axis and lies in xy plane)
        sensor_xyz(:,1) = receiver + [0;sensor_space/2;0]; % Add half sensor spacing to yp to give coords of L sensor
        sensor_xyz(:,2) = receiver + [0;-sensor_space/2;0]; % Subtract half sensor spacing from yp to give coords of R sensor
        % Check that both sensors are within room    
        if (min(min(sensor_xyz))<0)||(max(sensor_xyz(1,:))>=Lx)||(max(sensor_xyz(2,:))>=Ly)||(max(sensor_xyz(3,:))>=Lz)
            sensor_error='One or both sensors outside of room.';
            h=errordlg(sensor_error,terminate_msg);% Warn of impossible sensor location
            beep;
            uiwait(h); % Wait for user response
            sensor_outside=true;
        end;
        
    case {'mithrir','cipicir'}
        sensor_xyz(:,1) = receiver + [0;head_width/2;0]; % Add half ear spacing to yp to give coords of L sensor
        sensor_xyz(:,2) = receiver + [0;-head_width/2;0]; % Subtract half ear spacing from yp to give coords of R sensor
        % Check that both ears are within room    
        if (min(min(sensor_xyz))<0)||(max(sensor_xyz(1,:))>=Lx)||(max(sensor_xyz(2,:))>=Ly)||(max(sensor_xyz(3,:))>=Lz)
            sensor_error=strvcat('One or both sensors outside of room.');
            h=errordlg(sensor_error,terminate_msg);% Warn of impossible sensor location
            beep;
            uiwait(h); % Wait for user response
            sensor_outside=true;
        end;
end;
%--------------------- End of sensor location checks -----------------------------

B=realsqrt(1-A); % Calculate the frequency dependent pressure reflection coefficients (NB +ve for air to hard surface interface)
mean_B=mean(B,1); % Calculate the mean over frequency (ie. down the columns) of the reflection coefficients B

% Display the room geometry with receiver and source(s) locations as a 3D plot for visual confirmation
if alpha_F,
    alph=floor(64*mean_B'); %Opacity of each face proportional to average reflectivity NB TRANSPOSE
else, % Fix transparency for Ax1,Ax2,Ay1,Ay2,Az1,Az2
    alph=[32; 32; 32; 32; 64; 1]; % Opacity of faces fixed, opaque floor, clear ceiling, 50% transparent walls
end;
roomplot(room_size,source_xyz,receiver,sensor,sensor_xyz,head_width,alph,fig_loc);
%----------------------------- End of drawing the room --------------------------

if source_outside || sensor_outside, % Source or sensor outside room, so warning beep, and exit. 
    beep;
    return; % Force Exit
end
%-------------------------------- End of sanity checks ----------------------------------------

%-------------------- Data passed sanity checks so do the simulation -------------
fprintf(log_fid,'\n\n Speed of sound  c = %d m/s ',c); % Print c to the screen
%---------------- Set up the air absorption model ---------------------------
m_air=zeros(size(F_abs)); % Set air absorption to zero at octave frequencies
[RT60 MFP]= Norris_Eyring(c,room_size,A,F_abs,m_air); % Calculate a room reverberation time RT60 without air absorption
if air_F
    % The variation of air absorption coefficient "m" with humidity and frequency can be approximated by
    %    (from L E Kinsler et al "Fundamentals of acoustics", 3rd Ed., Wiley, New York, 1982.
    m_air = 5.5E-4*(50/humidity)*(F_abs/1000).^(1.7); % Valid Ranges are: relative humidity 20%< h <70%, 1500< F_abs <10000 Hz.
    % Calculate a room reverberation time RT60 to use in an estimate of a maximum impulse response length
    %(i.e. used to size H and set "order").
    [RT60_Air MFP_Air]= Norris_Eyring(c,room_size,A,F_abs,m_air); % Calculate a room reverberation time RT60 with air absorption
    if max(abs(RT60-RT60_Air)./RT60) < 0.01
        button = questdlg('Air absorption makes < 1% difference in this room, neglect?','Roomsim_run: Air Absorption Query','Yes','No','Yes');
        beep;
        if strcmp(button,'Yes')
            m_air=zeros(size(F_abs)); % Zero air absorption at standard octave frequencies
            air_F=0; % Clear the air flag (set false)
        else
            RT60=RT60_Air; % Use RT60 including air absorption
            MFP=MFP_Air; % Use Mean Free Path including air absorption
        end
    end
end
%---------------- End of air absorption set up ---------------
fprintf(log_fid,'\n\n Reverberation Times'); % Print to the screen
for ii=1:length(F_abs)
    fprintf(log_fid,'\n At frequency = %6.0f Hz RT60 = %8.4g s', F_abs(ii), RT60(ii)); % Print RT60 to the screen
end


%Estimate room break frequencies according to Everest p324 (4th Ed)
F1=c/max([Lx Ly Lz]); % Below F1 no resonant support for sound in room ie Lowest room mode
F2=11250*sqrt(0.0283*mean(RT60)/(Lx*Ly*Lz)); %Between F1 and F2 room modes dominate, between F2 and F3 diffraction and diffusion dominate
F3=4*F2; % Above F3 specular reflections and ray acoustics are valid

fprintf(log_fid,'\n\n Estimate of Room break frequencies (uses mean of above RT60''s),'); % Print to the screen
fprintf(log_fid,'\n F1 = %6.0f Hz. Lowest room mode, i.e. below F1 no resonant support for sound in room', F1);
fprintf(log_fid,'\n Between F1 and F2 room modes dominate'); % Print to the screen
fprintf(log_fid,'\n F2 = %6.0f Hz. Approximate cutoff (crossover) frequency', F2);
fprintf(log_fid,'\n Between F2 and F3 diffraction and diffusion dominate'); % Print to the screen
fprintf(log_fid,'\n F3 = %6.0f Hz. Above F3 specular reflections and ray acoustics are valid', F3);
%------------- Determine values for H_length and order ----------
if max(RT60)==0 % Detect anechoic case and set values to ensure full code path executed.
    order = [2; 2; 2]; 
    H_length = ceil(max(source_polar(1,:))*Fs_c); % Round up distance in samples between furthest source and receiver
    beep;
    msg_title='Roomsim_run: Anechoic room detected.';
    message=['Forcing order = ',num2str(max(order)),' and Impulse response length = ',num2str(H_length)];
    h=msgbox(message,msg_title,'warn');  %Warn & beep.
    beep;
    uiwait(h);% Wait for user to acknowledge
    
else
    if H_length < 0 % Estimate a reasonable maximum impulse response length
        H_length = ceil(max(RT60)*Fs); % H_length = longest reverberation time in samples (rounded up to integer)
    end
    
    if H_length > 15000 % Warn of long simulation time
        answer={};
        while isempty(answer),msg_title = 'Roomsim_run: Long impulse response';
            msg_title='Roomsim_run: Impulse response length shown may require a long simulation time';
            prompt = {'Enter a new impulse response length in samples (or keep current value):'};
            lines = 1;
            def = {num2str(H_length)}; %Default value shows current value
            answer = inputdlg(prompt,msg_title,lines,def,'on');
            beep;
        end; % Of while loop to disable inappropriate CANCEL button operation
        H_length=str2num(answer{1});
    end
    
    if order < 0 % then either Number of Reflections or impulse response length is used to limit order of reflections computed.
        if H_length < 0
            NoR = c*max(RT60)/MFP; % Estimate number of reflections
            order = ceil(NoR.^(1/3)); % Estimate order from Number of reflections
            order_x = order;
            order_y = order;
            order_z = order;
        else % Estimate order based on impulse response length (As Allen & Berkley)
            temp=c*H_length*T;
            NoR = temp/MFP; % Number of reflections
            order_x = ceil(temp/(2*Lx));
            order_y = ceil(temp/(2*Ly));
            order_z = ceil(temp/(2*Lz));
        end
        fprintf(log_fid,'\n\n Estimated number of reflections = %i', ceil(NoR)); % Print Number of reflections to the screen
    else % Accept user provided value of order
        order_x = order;
        order_y = order;
        order_z = order;    
    end
    
    order=[order_x; order_y; order_z]; % Pack up order into column vector
    
    if max(RT60)> 3 && max(order)>4 % Warn of long simulation time
        answer={};
        while isempty(answer),msg_title = 'Roomsim_run: Large value of Order';
            msg_title='Roomsim_run: Value of order may require a long simulation time';
            prompt = {'Enter a lower value for order (or keep current value) :'};
            lines = 1;
            def = {num2str(max(order))}; %Default value
            answer = inputdlg(prompt,msg_title,lines,def,'on');
            beep;
        end; % Of while loop to disable inappropriate CANCEL button operation
        temp=str2num(answer{1});
        order=min(order,temp); % Replace order values > temp with temp.
    end
    
    % Print order value to the screen
    fprintf(log_fid,'\n\n order_x = %i', order(1)); 
    fprintf(log_fid,'\n order_y = %i', order(2));
    fprintf(log_fid,'\n order_z = %i', order(3));
    
end
%------------------------- End of H_length and Order calculation -------------
fprintf(log_fid,'\n\n Sensor =  %s ',sensor); % Print receiver type
switch sensor % ---------- Identify receiver system---------------------
    case {'one_mic'}
        % Unpack sensor directionality and display to screen
        fprintf(log_fid,'\n min_azim_sensor = %6.2f deg',sensor_dir(1));
        fprintf(log_fid,'\n max_azim_sensor = %6.2f deg',sensor_dir(2));
        fprintf(log_fid,'\n min_elev_sensor = %6.2f deg',sensor_dir(3));
        fprintf(log_fid,'\n max_elev_sensor = %6.2f deg',sensor_dir(4));
    case {'two_mic'}
        % Unpack sensor directionality and display to screen
        fprintf(log_fid,'\n min_azim_sensor = [ %6.2f %6.2f ] deg',sensor_dir(1,1),sensor_dir(1,2));
        fprintf(log_fid,'\n max_azim_sensor = [ %6.2f %6.2f ] deg',sensor_dir(2,1),sensor_dir(2,2));
        fprintf(log_fid,'\n min_elev_sensor = [ %6.2f %6.2f ] deg',sensor_dir(3,1),sensor_dir(3,2));
        fprintf(log_fid,'\n max_elev_sensor = [ %6.2f %6.2f ] deg',sensor_dir(4,1),sensor_dir(4,2));
end;

%--------- Display the reverberation time as a plot of RT60 against frequency ------------
roomplot_RT60(F_abs, RT60,fig_loc);
drawnow; %Force completion of all previous figure drawing before continuing
%--------------------------------------------------------------------------------
 
%Run the frequency dependent absorption model of the room with the parameter values as set above.
[source,p_isource,H,HRTF,n_sources] = roomsim_core(c,humidity,Fs,room_size,source_xyz,receiver...
    ,sensor,sensor_xyz,sensor_dir,S_No,F_abs,B,air_F,smooth_F,Fc_HP,dist_F,order,H_length,H_filename);

%----------------------- Plot the impulse responses to each sensor ---------------
roomplot_imp(Fs,H,sensor,fig_loc);

%----------------------- Plot the frequency responses to each sensor ---------------
roomplot_magf(Fs,H,HRTF,sensor,fig_loc);

%----------------------- 2D Plot of the room and the image rooms---------------
if plot_F2 % Display the receiver, source, image sources and image rooms as a 2D plan.
    n_images=size(p_isource,1); % Number of image sources
%     if n_images >300,
%         wait_msg = ['Processing 2D Plot, please wait...'];% Keep the user informed about progress
% %         msg_title='Roomplot_2D';
%         h_plot2D= waitbar(0.2,wait_msg); %Display wait message
%     end;
    roomplot_2D(c,Fs,room_size,source_xyz,receiver,sensor,sensor_xyz,sensor_dir...
        ,H_length,source,p_isource,fig_loc);
%     if n_images >300 & exist('h_plot2D','var'),
%         close(h_plot2D);
%     end; % Close the wait message box if it was opened
end;

%----------------------- 3D Plot of the room and the image sources---------------
if plot_F3 % Display the receiver, source, image sources and image rooms in 3D
    n_images=size(p_isource,1); % Number of image sources
%     if n_images >100,
%         wait_msg = ['Processing 3D Plot, please wait...'];% Keep the user informed about progress
% %         msg_title='Roomplot_3D';
%         h_plot3D=waitbar(0.2,wait_msg); %Update progress
%     end;
    roomplot_3D(room_size,source_xyz,receiver,sensor,sensor_xyz,sensor_dir,source,p_isource,fig_loc);
%     if n_images >100 & exist('h_plot3D','var'),
%         close(h_plot3D);
%     end; % Close the wait message box if it was opened
end;

%----------------- Save data for plotting to a MATLAB loadable file-----------------
drawnow; %Force completion of all figure drawing before continuing
beep;
file_spec='PLOT_Roomsim.mat';
title='Roomsim_run: Save data as a MAT file for later plotting';
dot_ext=[];
while isempty(dot_ext),
    [name path] = uiputfile(file_spec, title); %Display the dialogue box
    if ~any(name), % File select was cancelled, exit the while loop.
        fprintf(log_fid,'\n\n Roomsim_run: Plot data has not been saved'); % Print to log file (avoids wait for user response)
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
                save([path name], 'c','Fs','room_size','source_xyz','receiver','sensor','sensor_xyz','sensor_dir','H_length','source','p_isource','H','HRTF'); % Save relevant plot variables to a ".mat" file
                fprintf(log_fid,'\n\n Roomsim_run: Plot data has been saved to %s', [path name]); % Print to log file (avoids wait for user response)
            otherwise,
                h=errordlg('Data file extension not recognised, exiting','savef error');
                beep;
                uiwait(h);
                return;
        end;
    end;
end; % Of while loop to trap missing .mat extension

%------------------------------ End of roomsim_run.m ------------------------------
