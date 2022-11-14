function [source,p_isource,H,HRTF,n_sources]=roomsim_core(c,humidity,Fs,room_size,source_xyz,receiver,sensor,sensor_xyz,sensor_dir,S_No...
    ,F_abs,B,air_F,smooth_F,Fc_HP,dist_F,order,H_length,H_filename)
%Usage: [source,p_isource,H,HRTF,n_sources]=roomsim_core(c,humidity,Fs,room_size,source_xyz,receiver,sensor,sensor_xyz,sensor_dir,S_No...
%     ,F_abs,B,air_F,smooth_F,Fc_HP,dist_F,order,H_length,H_filename);
% roomsim_core.m is a MATLAB m-code implementation of a mirror image model of the impulse response from an 
%   omni-directional source to a receiver (sensor, two sensors or two ears of a head) in a "shoebox" room.
% It incorporates frequency dependent absorption at the reflective surfaces and in the airspace of the room.
% The simulation of the head utilises Head Related Transfer Function (HRTF) data
%   actually Head Related Impulse Response (HRIR) data provided from measurements made either;
%   1)on a Kemar mannequin at MIT. (http://xenia.media.mit.edu/~kdm/hrtf.html ,Dec 2002)
%   or
%   2)on real subjects and a Kemar at University of California, Davis. (http://interface.cipic.ucdavis.edu ,Dec 2002).
% The kernel of the image method in Roomsim is derived from the Fortran program reported by Allen and Berkley,1979 [A6]).
% -------------------------------------------------------------------------
% Copyright (C) 2003  Douglas R Campbell and Kalle Palomaki
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
%-------------------------------------------------------------
% The original version of this MATLAB implementation was written by Kalle Palomaki 
%   (Helsinki University of Technology, Finland), while working (with partial support 
%   from a Finnish Tekniikan edistamissaatio grant) on the EC TMR SPHEAR project at the University of Sheffield.
%   The code was subsequently modified and incorporated into this menu driven program by Douglas R. Campbell 
%   (University of Paisley, Scotland, UK) during Sept2002/Feb 2003.
%--------------------------------------------------------------------------------------------------------------
% References for image method program code:
% (A5) P M Peterson, "Simulating the response of multiple microphones to a single acoustic source in a reverberant room", JASA 80(5),1986, 1527-1529.
% (A6) J B Allen and D A Berkley, "Image method for efficiently simulating small-room acoustics", JASA 65(4), 1979, 943-950.
%-----------------------------------------------------------------------------------------------------------------
% Functions called: getNearestMITpulse.m, getNearestUCDpulse.m
global log_fid; % Identifier of logfile
global hrir_l hrir_r;% Globals for CIPIC impulse responses left and right
global H3D; % Global for MIT impulse response data

%Initialise constants
degtorad=pi/180; % Conversion factor degrees to radians
radtodeg=180/pi; % Conversion factor radians to degrees
Two_pi=2*pi; % Compute here for efficiency
T=1/Fs; % Sampling Period
nyquist=Fs/2; % Half sampling frequency
Fs_c = Fs/c; % Samples per metre
F_abs=[0 F_abs nyquist]; % Extend F_abs to include the 0 and Fs/2 Hz values

% Unpack room dimensions
Lx=room_size(1); % Length
Ly=room_size(2); % Width
Lz=room_size(3); % Height

%Unpack source(s) coordinates
x=source_xyz(1,:);
y=source_xyz(2,:);
z=source_xyz(3,:);

%Unpack surface reflectance
bx1=B(:,1);
bx2=B(:,2);
by1=B(:,3);
by2=B(:,4);
bz1=B(:,5);
bz2=B(:,6);

%Unpack order
order_x = order(1);
order_y = order(2);
order_z = order(3);

%Unpack receiver reference coordinates (Head or sensor(s))
xp=receiver(1);
yp=receiver(2);
zp=receiver(3);

%------------- Specific setups for each sensor case -------------------------------
switch sensor % ---------- Identify receiver system---------------------
    case {'one_mic','two_mic'}
        % Unpack sensor directionality and display to screen
        min_azim_sensor=sensor_dir(1,:);
        max_azim_sensor=sensor_dir(2,:);
        min_elev_sensor=sensor_dir(3,:);
        max_elev_sensor=sensor_dir(4,:);
end;

%------------------------- Set up Air attenuation model ------------
if air_F  % Include the absorption due to air
    % Estimate the frequency dependent pressure absorption coeff for air, m= 2*alpha (neper/m)
    m_air = 5.5E-4*(50/humidity)*(F_abs/1000).^(1.7);
    atten_air= exp(-0.5*m_air); %Compute the attenuation factor for one metre travelled in air
end;

%------------------------- Set up Smoothing Filter ------------
% For use when HRTF is not present (Ref A5) to take account of impulses not centered on a sampling instant
Fc = 0.9*nyquist; %Smoothing Filter Cutoff marginally below Fs/2, >40dB attenuation close to Fs/2.
N_smooth = 32; % Order of FIR smoothing filter
Fc_Fs=Fc.*T; % Compute here for efficiency
Two_Fc=2.*Fc; % Compute here for efficiency
Tw = N_smooth*T; % Window duration (seconds)
Tw_inv=1/Tw; % Invert here for efficiency
Two_pi_Tw=Two_pi.*Tw_inv; % Compute here for efficiency
t=[-Tw/2:T:Tw/2]'; % Filter time window NB column vector symmetrical about t=0

%------------------------- Set up High-pass Filter ------------
%Used in Allen & Berkley (and probably CoolEdit).This appears to be an attempt to overcome an 
%accumulating DC offset due (I guess) to absence of pressure equalisation (i.e. leakage) in the "ideal" room.
%Probably not necessary if the sampling frequency is sufficently high e.g. Fs=10*fh
if Fc_HP > 0,
    [b_HP,a_HP] = hi_pass2(Fc_HP,Fs); % Compute coefficients of the Allen & Berkley 2nd order high-pass filter
    %   load HP_coeffs; IF active this code uses the 4th order Butterworth High-pass filter of the utilities suite.
    %   b_HP=B;
    %   a_HP=A;
end;
%------------------------------ General Set ups for the image model ------------------
% isource_ident codes the eight permutations of x+/-xp,y+/-yp,z+/-zp (the source to receiver vector components)
%   where [-1 -1 -1] identifies the parent source.
isource_ident=[-1 -1 -1; -1 -1 1; -1 1 -1; -1 1 1; 1 -1 -1; 1 -1 1; 1 1 -1; 1 1 1];
surface_coeff=[0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1]; % Includes/excludes bx,by,bz depending on 0/1 state.

%Allocate power indice modifiers here for efficiency
qq=surface_coeff(:,1); %  for bx1
jj=surface_coeff(:,2); %  for by1
kk=surface_coeff(:,3); %  for bz1

%Dimension and clear arrays
n_sources = size(source_xyz,2); % Number of parent sources (size finds number of columns ie sources)
n_isources = (2*order_x+1)*(2*order_y+1)*(2*order_z+1)*8; %Maximum number of image sources
p_isource = zeros(n_isources,n_sources,2); % Allocate array for image source maximum strength value, one per sensor (DIM(3)=1 is L.
H = zeros(H_length, 2, n_sources); % Allocate array for L&R impulse responses from each source
HRTF=[]; %Declare array HRTF
source = zeros(3,n_isources,n_sources); % Allocate array for image source co-ordinates

F_abs=F_abs./nyquist; %Normalise the standard absorption frequency range for surfaces, (0 to 1) = (0 to Fs/2)

% Set up parameters for converting surface frequency response to impulse response
N_refl=128; %Required length of FIR filter modelling impulse response of surface(+air)
Half_I=fix(N_refl./2); % Shift due to FIR filter model
interp_x=[min(F_abs): (max(F_abs)-min(F_abs))./Half_I :max(F_abs)]; %Generate the interpolation grid
window = 0.5.*(1 - cos(2*pi*(0:N_refl).'./N_refl));%Compute the Hann window

%-------------------- Calculate the impulse responses from each primary source location ------------------------------
msg_title='Roomsim_core:';
wait_msg_img = ['  Generating Images of ' num2str(n_sources) ' Primary Source(s), please wait...'];
wait_msg_imp = '   Computing Impulse Responses, Please wait...';
h_parent = waitbar(0.1,wait_msg_img,'name',msg_title); 

for ps=1:n_sources %Compute for each parent source
    % Keep the user informed about progress. Increment the waitbar a little to show its started
    beep;

    temp_source = zeros(3,n_isources); % Allocate array for temporary image source co-ordinates
    dist = zeros(1,n_isources); % Allocate array for distances from image source to receiver
    refl = zeros(6,n_isources); % Allocate array for surface reflection impulse amplitude
    xx=zeros(8,1); yy=zeros(8,1); zz=zeros(8,1); %Allocate work arrays
    hrir=[]; % Declare HRTF impulse response array
    %----------Compute these here for efficiency ------------------------
    xx=isource_ident(:,1)*x(ps); %partial x coord of image.
    yy=isource_ident(:,2)*y(ps); %partial y coord of image.
    zz=isource_ident(:,3)*z(ps); %partial z coord of image. 
    xx_yy_zz=[xx yy zz]';
    Two_Lx=2*Lx; Two_Ly=2*Ly; Two_Lz=2*Lz;
    
    %------------- Calculate the amplitude and location of each of the image sources ------------------   
    % Computes the frequency dependent surface reflection and coordinates and distance for each image
    % Overall surface reflection in each octave band (1:6), coordinates, and distance for each image source n_images
    % computed using refl(:,n_images)=bx1(:)^abs(n-q)*bx2(:)^abs(n)*by1(:)^abs(l-j)*by2(:)^abs(l)*bz1(:)^abs(m-k)*bz2(:)^abs(m)
    % with partials of this expression pre-computed for efficiency.
    
    waitbar(ps/n_sources,h_parent); % Increment the waitbar
%     tic;  % Start TIMER
    
    n_images=0; %Clear n_images, used to count the number of significant images of each parent source
    
    for n=-order_x:order_x,
        bx2_abs_n=bx2.^abs(n); %Compute here for efficiency
        Two_n_Lx=n*Two_Lx; %Compute here for efficiency
        
        for l=-order_y:order_y,
            bx2y2_abs_nl=bx2_abs_n.*(by2.^abs(l)); %Compute here for efficiency
            Two_l_Ly=l*Two_Ly; %Compute here for efficiency
            
            for m=-order_z:order_z,
                bx2y2z2_abs_nlm=bx2y2_abs_nl.*(bz2.^abs(m)); %Compute here for efficiency
                Two_m_Lz=m*Two_Lz; %Compute here for efficiency
                
                Two_nlm_Lxyz = [Two_n_Lx; Two_l_Ly; Two_m_Lz]; %Concatenate here for efficiency
                
                for permu=1:8,
                    n_images=n_images+1; %Accumulate count of the image sources
                    % calculate xyz coordinates of image source n_images of parent source ps
                    temp_source(:,n_images)=Two_nlm_Lxyz - xx_yy_zz(:,permu);
                    % Calculate distance (using Pythagoras in 3D) from receiver to each parent and image source.
                    dist(n_images) = norm((temp_source(:,n_images)-receiver(:)),2); %Distance (m) between source and receiver reference point                       
                    % Calculate delay in samples between receiver reference point and each image source
                    delay = Fs_c*dist(n_images); % delay in samples = (Samples per metre)*Distance
                    if delay <= H_length,
                        refl(:,n_images)=bx1.^abs(n-qq(permu)).*by1.^abs(l-jj(permu)).*bz1.^abs(m-kk(permu)).*bx2y2z2_abs_nlm;                    
                        if max(refl(:,n_images)) < 1E-3, % (NB refl always +ve for air to hard surface, otherwise need abs here)
                            n_images=n_images-1; % Delete image sources with a max reflection coeff below 1*10^-3
                        end
                    else
                        n_images=n_images-1; % Delete image sources with a delay > H_length
                    end; % of decimation of low level or distant sources
                end; % of permu counter loop
                
            end; % of m counter loop
        end; %of l counter loop
    end; % of n counter loop and generation of n_images image source(s)
    
    %--------------------------- Compute the complete impulse response, Primary source(s) to Receiver ---------------------   
    
    temp_source = temp_source(:,1:n_images); % Re-Allocate array for image source co-ordinates (discard trailing zero values)
    dist = dist(1:n_images); % Re-Allocate array for distances from image source to receiver (discard trailing zero values)
    refl = refl(:,1:n_images); % Re-Allocate array for surface reflection impulse amplitude (discard trailing zero values)
    azim = zeros(1,n_images); %Allocate array for image source azimuths
    elev = zeros(1,n_images); %Allocate array for image source elevations
    b_refl = zeros(n_images,8); %Allocate array for (partial) surface reflectancies
    Hp = zeros(H_length,2); %Allocate array Hp (two column stereo format for sound files), make it long enough to hold impulse response length H_length.
    
    inv_dist=ones(size(dist)); % Allocate and initialise array for 1/dist. Compute here for efficiency.
    if dist_F, % Apply distance (ie. 1/R) attenuation. Allows easy removal of effect of distance for checking.
        % By setting dist <= 1m to 1m an attenuation of 0 dB is applied to sources within 1 metre of receiver
        % Avoids amplification when 0<dist<1, and divide by zero when dist=0 i.e. source coincident with receiver
        inv_dist(dist>1) = 1./dist(dist>1); % Using logical(dist>1) for array indexing
    end;
    
    refl = [refl(1,:); refl(:,:); refl(6,:)]; % Pad refl to same length as F_abs by estimating the 0 and Fs/2 Hz values
    for ind=1:8,
        b_refl(:,ind) = (refl(ind,:).*inv_dist).';% Scale to take account of image source distance NB Transposing here
    end;%         
    
    % Keep the user informed about progress. Increment the waitbar a little to show its started
    h_wbar = waitbar(0.1,wait_msg_imp,'name',msg_title);
    inv_n_images = 1/n_images; % Compute here for efficiency
    step = ceil(n_images/5); % Waitbar step size, 5 steps displayed
    
    for is = 1:n_images % for each of the n_images image sources
        
        if mod(is,step)==0, % Update the waitbar in steps of size step
            waitbar(is*inv_n_images,h_wbar); % Increment the waitbar
        end;
        
        if air_F  % Include the absorption due to air
            %Compute a set of distance dependent pressure reflection coefficients for air
            b_refl(is,:)=b_refl(is,:).*(atten_air.^dist(is)); % NB. Assumes distance to reference point is appropriate even for two sensors
        end
        
        %--------------- Transform from frequency response to impulse response -------------------------------
        % Replacement for fir2 to avoid requiring users to install SIGNAL toolbox (also faster)
        yy1 = interp1(F_abs',b_refl(is,:),interp_x); % Interpolate to estimate the values at the dense grid points
        yy3 = [yy1 conj(fliplr(yy1(1:(Half_I))))]; % Real data has symmetric Mag and anti-sym phase about Nyquist
        yy3 = yy3(1:length(yy3)-1)'; % Transpose to column and discard last data point to make periodic spectrum  
        total_refl = real(ifft(yy3,length(yy3))); % IFFT to calculate impulse
        total_refl = window.*[total_refl(Half_I+1:N_refl); total_refl(1:Half_I+1)]; % Make the impulse realisable (half length shift) and Hann window it 
        %------------------------------------------------------------------------------------------------------
        
        % Ignore sources with impulse response peak magnitudes below -60dB (1/1000) unless anechoic case  
        if (n_images==1) | (max(abs(total_refl)) > 1E-3)
            switch sensor % ---------- Select receiver system---------------------
                % Calculate the Time of Arrival (ToA) for the signal from each image source and accumulate the impulse response
                case 'one_mic'
                    dist_m = dist(is); % norm((temp_source(:,is)-sensor_xyz),2); %Distance (m) between source and single sensor (receiver)
                    delay = Fs_c*dist_m; % delay in samples = (Samples per metre)*Distance
                    rdelay=round(delay); % Extract integer delay (concatenated later with impulse response)
                    if dist_m < eps, % Source coincident with sensor
                        elev = 0; % angles irrelevant
                        azim = 0;
                    else, % angles important
                        elev = radtodeg*asin((temp_source(3,is)-zp)/dist_m); % Calculate elevation -90<= elev <=90 deg
                        azim = radtodeg*angle((temp_source(1,is)- xp)+j*(temp_source(2,is)- yp)); %Calculate -180 >= azimuth >= 180 deg
                    end;
                    if (azim >= min_azim_sensor)&&(azim <=max_azim_sensor)&&(elev >= min_elev_sensor)&&(elev <=max_elev_sensor) %  Check elevation & azimuth visibility
                        if smooth_F % Include smoothing filter
                            t_Td=t-(delay-rdelay).*T; % Take account of fractional delay  -0.5 < D < +0.5 sample period
                            hsf=Fc_Fs*(1+cos(Two_pi_Tw.*(t_Td))).*sinc(Two_Fc*(t_Td)); % Compute delayed filter impulse response for Left sensor
                            %----- Do the convolutions using the FILTER primitive for speed, shortest sequence is first parameter.
                            total_refl=[total_refl; zeros(size(hsf,1)-1,1)]; % Append filter_length-1 zeros so convolution length is data_length+impulse_length-1
                            h = filter(hsf, 1, total_refl); % Convolve Left channel signals
                        else % No smoothing
                            h = total_refl;
                        end
                        
                        p_isource(is,ps,:)=max(abs(h)); % Max of impulse response from this image (used for plotting image source strength)
                        
                        len_h=length(h); %
                        adjust_delay = rdelay - ceil(len_h/2); % Half length shift to remove delay due to impulse response
                        
                        %--- Accumulate the impulse responses from each image source within an array of length H_length ---
                        if adjust_delay < 0
                            start_index_Hp = max(adjust_delay+1,1);
                            stop_index_Hp = min(adjust_delay+1+len_h,H_length);
                            start_index_h = max(-adjust_delay,1);
                            stop_index_h = min(H_length-adjust_delay,len_h);
                        else
                            start_index_Hp = max(adjust_delay+2,1);
                            stop_index_Hp = min(adjust_delay+len_h+1,H_length);
                            start_index_h = max(-adjust_delay,1);
                            stop_index_h = min(H_length-adjust_delay-1,len_h);
                        end;
                        Hp(start_index_Hp:stop_index_Hp,1)= Hp(start_index_Hp:stop_index_Hp,1) + h(start_index_h:stop_index_h); %Add whole or part of impulse response
                        
                        Hp(:,2)=Hp(:,1); % Only one sensor so copy the impulse response to right channel
                    end; % of Check visibility etc
                    
                case 'two_mic'
                    %Microphone pair centered on (xp,yp,zp) lying in xy plane and aligned with x axis separated by a spacing of sensor_space metres
                    dist_m(1) = norm(temp_source(:,is)-sensor_xyz(:,1),2); %Distance (m) between source and receiver left sensor
                    dist_m(2) = norm(temp_source(:,is)-sensor_xyz(:,2),2); %Distance (m) between source and receiver right sensor
                    delay = Fs_c.*dist_m; % L & R delays in samples = (Samples per metre)*Distance
                    rdelay=round(delay); % Extract L & R integer delays (concatenated later with impulse response)
                    
                    elev = radtodeg.*asin((temp_source(3,is)-sensor_xyz(3,:))./dist_m); % Calculate image source elevations from L&R sensors -90<= elev <=90 deg
                    azim = radtodeg.*angle((temp_source(1,is)- sensor_xyz(1,:))+j.*(temp_source(2,is)- sensor_xyz(2,:))); %Calculate -180 >= L&R azimuths >= 180 deg
                    
                    % ---------------------------- Left sensor ------------------------------
                    if (azim(1) >= min_azim_sensor(1))&&(azim(1) <=max_azim_sensor(1))&&(elev(1) >= min_elev_sensor(1))&&(elev(1) <=max_elev_sensor(1))%  Check elevation & azimuth visibility
                        %ignore sources not visible to L sensor ie. min_azim_sensor > azim > max_azim_sensor deg, min_elev_sensor > elev > max_elev_sensor deg
                        if smooth_F % Include smoothing filter
                            t_Td=t-(delay(1)-rdelay(1)).*T; %Take account of fractional delay  -0.5 < D < +0.5 sample period
                            hsf=Fc_Fs.*(1+cos(Two_pi_Tw.*(t_Td))).*sinc(Two_Fc.*(t_Td)); % Compute delayed filter impulse response for Left sensor
                            %----- Do the convolutions using the FILTER primitive for speed.
                            total_refl_L=[total_refl; zeros(length(hsf)-1,1)]; % Append filter_length-1 zeros so convolution length is data_length+impulse_length-1
                            h_L = filter(hsf, 1, total_refl_L); % Convolve Left channel signals
                        else % No smoothing
                            h_L = total_refl;
                        end
                        p_isource(is,ps,1)=max(abs(h_L)); % Max of impulse response in L channel (used for plotting image source strength)
                        
                        len_h_L=length(h_L); %
                        adjust_delay_L = rdelay(1) - ceil(len_h_L/2); % Half length shift of impulse response
                        
                        %---  Accumulate the L impulse responses from each image source within an array of length H_length ---
                        
                        if adjust_delay_L < 0
                            start_index_Hp = max(adjust_delay_L+1,1);
                            stop_index_Hp = min(adjust_delay_L+1+len_h_L,H_length);
                            start_index_h = max(-adjust_delay_L,1);
                            stop_index_h = min(H_length-adjust_delay_L,len_h_L);
                        else
                            start_index_Hp = max(adjust_delay_L+2,1);
                            stop_index_Hp = min(adjust_delay_L+len_h_L+1,H_length);
                            start_index_h = max(-adjust_delay_L,1);
                            stop_index_h = min(H_length-adjust_delay_L-1,len_h_L);
                        end;
                        Hp(start_index_Hp:stop_index_Hp,1)= Hp(start_index_Hp:stop_index_Hp,1) + h_L(start_index_h:stop_index_h); %Add whole or part of impulse response
                    end; % of Check visibility Left etc
                    
                    % ---------------------------- Right sensor -------------------------------
                    if (azim(2) >= min_azim_sensor(2))&&(azim(2) <=max_azim_sensor(2))&&(elev(2) >= min_elev_sensor(2))&&(elev(2) <=max_elev_sensor(2))%  Check elevation & azimuth visibility
                        %ignore sources not visible to R sensor ie. min_azim_sensor > azim > max_azim_sensor deg, min_elev_sensor > elev > max_elev_sensor deg
                        if smooth_F % Include smoothing filter
                            t_Td=t-(delay(2)-rdelay(2)).*T; %Take account of fractional delay  -0.5 < D < +0.5 sample period
                            hsf=Fc_Fs.*(1+cos(Two_pi_Tw.*(t_Td))).*sinc(Two_Fc.*(t_Td)); % Compute delayed filter impulse response for Right sensor
                            %----- Do the convolutions using the FILTER primitive for speed.
                            total_refl_R=[total_refl; zeros(length(hsf)-1,1)]; % Append filter_length-1 zeros so convolution length is data_length+impulse_length-1
                            h_R = filter(hsf, 1, total_refl_R); % Convolve Right channel signals
                        else % No smoothing
                            h_R = total_refl;
                        end
                        p_isource(is,ps,2)=max(abs(h_R)); % Max of impulse response in R channel
                        
                        len_h_R=length(h_R); %
                        adjust_delay_R = rdelay(2) - ceil(len_h_R/2); % Half length shift of impulse response
                        
                        %---  Accumulate the R impulse responses from each image source within an array of length H_length ---
                                                
                        if adjust_delay_R < 0
                            start_index_Hp = max(adjust_delay_R+1,1);
                            stop_index_Hp = min(adjust_delay_R+1+len_h_R,H_length);
                            start_index_h = max(-adjust_delay_R,1);
                            stop_index_h = min(H_length-adjust_delay_R,len_h_R);
                        else
                            start_index_Hp = max(adjust_delay_R+2,1);
                            stop_index_Hp = min(adjust_delay_R+len_h_R+1,H_length);
                            start_index_h = max(-adjust_delay_R,1);
                            stop_index_h = min(H_length-adjust_delay_R-1,len_h_R);
                        end;
                        Hp(start_index_Hp:stop_index_Hp,2)= Hp(start_index_Hp:stop_index_Hp,2) + h_R(start_index_h:stop_index_h); %Add whole or part of impulse response
                    end; % of Check visibility Right etc
                    
                case 'mithrir'
                    delay = Fs_c*dist(is); % delay in samples = (Samples per metre)*Distance
                    rdelay=round(delay); % Extract integer delay (concatenated later with impulse response)
                    
                    %MIT Kemar is present so calculate elevation and azimuth of each image source
                    elev = radtodeg*asin((temp_source(3,is)-zp)/dist(is)); % Calculate elevation  -90<= elev <=90 deg             
                    if (elev >= -40)&(elev <=90) % Check visibility, ignore sources 220 <elev < 320 deg and elev > 90 deg
                        len_hh=length(total_refl); %
                        adjust_delay = rdelay - ceil(len_hh/2); % Half length shift of impulse response
                        
                        azim = radtodeg*angle((temp_source(1,is)- xp)+j*(temp_source(2,is)- yp)); %Calculate -180 >= azimuth >= 180 deg

                        % Quantise the elevations and azimuths to nearest angle in the MIT set (-40 >= elevation >= +90 deg, -180 >= azimuth >= +180 deg)
                        %   and get the HRIR for this i'th image source NB MIT HRIR starts with an ~3 sample offset at elev=0deg azim=+/- 90deg
                        [hrir, elerr, azerr] = getNearestMITpulse(elev,azim,H3D); % Get the HRIR for this i'th image source
                        
                        % Convolve the reflection and head related impulse responses (hrir) using the FILTER primitive for speed.
                        hrir_ext=[hrir; zeros(length(total_refl)-1,2)]; % Append filter_length-1 zeros so convolution length is data_length+impulse_length-1
                        h(:,1)=filter(total_refl,1,hrir_ext(:,1)); % Left ear
                        h(:,2)=filter(total_refl,1,hrir_ext(:,2)); % Right ear                        
                        p_isource(is,ps,:)=max(abs(h),[],1); % Max's of impulse response in L & R channels (used for plotting image source strength)
                        
                        len_h=length(h(:,1)); %
                        %--- Accumulate the L impulse responses from each image source within an array of length H_length ---
                        if adjust_delay < 0
                            start_index_Hp = max(adjust_delay+1,1);
                            stop_index_Hp = min(adjust_delay+1+len_h,H_length);
                            start_index_h = max(-adjust_delay,1);
                            stop_index_h = min(H_length-adjust_delay,len_h);
                        else
                            start_index_Hp = max(adjust_delay+2,1);
                            stop_index_Hp = min(adjust_delay+len_h+1,H_length);
                            start_index_h = max(-adjust_delay,1);
                            stop_index_h = min(H_length-adjust_delay-1,len_h);
                        end;
                        Hp(start_index_Hp:stop_index_Hp,1)= Hp(start_index_Hp:stop_index_Hp,1) + h(start_index_h:stop_index_h,1); %Add whole or part of impulse response
                        
                        len_h=length(h(:,2)); %
                        %--- Accumulate the R impulse responses from each image source within an array of length H_length ---
                        if adjust_delay < 0
                            start_index_Hp = max(adjust_delay+1,1);
                            stop_index_Hp = min(adjust_delay+1+len_h,H_length);
                            start_index_h = max(-adjust_delay,1);
                            stop_index_h = min(H_length-adjust_delay,len_h);
                        else
                            start_index_Hp = max(adjust_delay+2,1);
                            stop_index_Hp = min(adjust_delay+len_h+1,H_length);
                            start_index_h = max(-adjust_delay,1);
                            stop_index_h = min(H_length-adjust_delay-1,len_h);
                        end;
                        Hp(start_index_Hp:stop_index_Hp,2)= Hp(start_index_Hp:stop_index_Hp,2) + h(start_index_h:stop_index_h,2); %Add whole or part of impulse response
                        
                        if temp_source(:,is)==source_xyz(:,ps), % image is co-incident with its parent source
                            HRTF(:,:,ps)=hrir(:,:); %Save hrir's for each parent source direction in HRTF for plotting
                        end;
                        
                    end; % of Check visibility etc
                    
                case 'cipicir'
                    delay = Fs_c*dist(is); % delay in samples = (Samples per metre)*Distance
                    rdelay=round(delay); % Extract integer delay (concatenated later with impulse response)
                    
                    %CIPIC Head is present so calculate elevation and azimuth of each image source
                    elev = radtodeg*asin((temp_source(3,is)-zp)/dist(is)); % Calculate elevation for visibility test -90<= elev <=90 deg                
                    if (elev >= -45)|(elev <=-129.375) % Check visibility, ignore sources 230.625 <elev < 315 deg
                        len_hh=length(total_refl); %
                        adjust_delay = rdelay - ceil(len_hh/2); % Half length shift of impulse response
                        
                        elev = radtodeg*angle((temp_source(1,is)- xp)+j*(temp_source(3,is)- zp)); %Calculate -180 >= elevation >= +180 deg
                        azim = radtodeg*asin((temp_source(2,is)-yp)/dist(is)); % Calculate -90 >= azimuth >= +90 deg                   
                        % Quantise the elevations and azimuths to nearest angle in the CIPIC set (-80 >= elevation >= +80 deg, -45 >= azimuth >= +230 deg)
                        %   and get the HRIR for this i'th image source NB CIPIC HRIR starts with an ~30 sample offset
                        [hrir(:,1), azerr, elerr] = getNearestUCDpulse(-azim,elev,hrir_l);% NB azimuth for CIPIC measurements was +ve CW so invert sign.
                        [hrir(:,2), azerr, elerr] = getNearestUCDpulse(-azim,elev,hrir_r);% NB azimuth for CIPIC measurements was +ve CW so invert sign.
                        
                        % Convolve the reflection and head related impulse responses (hrir) using the FILTER primitive for speed.
                        hrir_ext=[hrir; zeros(length(total_refl)-1,2)]; % Append filter_length-1 zeros so convolution length is data_length+impulse_length-1
                        h(:,1)=filter(total_refl,1,hrir_ext(:,1)); % Left ear
                        h(:,2)=filter(total_refl,1,hrir_ext(:,2)); % Right ear
                        
                        p_isource(is,ps,:)=max(abs(h),[],1); % Max's of impulse response in L & R channels (used for plotting image source strength)
                        
                        len_h=length(h(:,1)); %
                        %--- Accumulate the L impulse responses from each image source within an array of length H_length ---
                        if adjust_delay < 0
                            start_index_Hp = max(adjust_delay+1,1);
                            stop_index_Hp = min(adjust_delay+1+len_h,H_length);
                            start_index_h = max(-adjust_delay,1);
                            stop_index_h = min(H_length-adjust_delay,len_h);
                        else
                            start_index_Hp = max(adjust_delay+2,1);
                            stop_index_Hp = min(adjust_delay+len_h+1,H_length);
                            start_index_h = max(-adjust_delay,1);
                            stop_index_h = min(H_length-adjust_delay-1,len_h);
                        end;
                        Hp(start_index_Hp:stop_index_Hp,1)= Hp(start_index_Hp:stop_index_Hp,1) + h(start_index_h:stop_index_h,1); %Add whole or part of impulse response
                        
                        len_h=length(h(:,2)); %
                        %--- Accumulate the R impulse responses from each image source within an array of length H_length ---
                        if adjust_delay < 0
                            start_index_Hp = max(adjust_delay+1,1);
                            stop_index_Hp = min(adjust_delay+1+len_h,H_length);
                            start_index_h = max(-adjust_delay,1);
                            stop_index_h = min(H_length-adjust_delay,len_h);
                        else
                            start_index_Hp = max(adjust_delay+2,1);
                            stop_index_Hp = min(adjust_delay+len_h+1,H_length);
                            start_index_h = max(-adjust_delay,1);
                            stop_index_h = min(H_length-adjust_delay-1,len_h);
                        end;
                        Hp(start_index_Hp:stop_index_Hp,2)= Hp(start_index_Hp:stop_index_Hp,2) + h(start_index_h:stop_index_h,2); %Add whole or part of impulse response
                        
                        if temp_source(:,is)==source_xyz(:,ps),
                            HRTF(:,:,ps)=hrir(:,:); %Save hrir's for each parent source direction in HRTF for plotting
                        end;
                    end; % of Check visibility etc
                    
                otherwise
                    h_err=errordlg('Unknown receiver set up. Program will terminate.','Roomsim Error');
                    beep;
                    uiwait(h_err);
                    return
                    
            end;%-------------------- End of Receiver system (Switch Case) select ------------------------------- 
        end;%of trap low amplitude image sources
    end;% of is counter loop for number of image sources
    
    close(h_wbar); % Close the current waitbar displaying progress through the image sources
%     toc;   % Stop TIMER
    
    fprintf(log_fid,'\n\n Number of Images = %i ',n_images); % Print to the logfile
    
    len_source(ps) = size(temp_source,2); % Keep track of length of the impulse response for each primary source
    source(:,1:len_source(ps),ps) = temp_source(:,:); % Save image source co-ordinates in array source for plotting
    
    if Fc_HP > 0,
        Hp(:,1) = filter(b_HP,a_HP,Hp(:,1)); % ARMA High-pass filter the left channel impulse response (Allen and Berkley used Fc_HP=100 Hz)
        Hp(:,2) = filter(b_HP,a_HP,Hp(:,2)); % Do the right channel also.
        %  Hp = Hp-filter(ones(1,64)/64,1,Hp); % MA High-pass filter possible alternative if delay needs compensated
    end;
    
    H(:,:,ps)=Hp(:,:); % Save accumulated impulse response due to each parent source in H for plotting 
    
end;% of ps counter loop for number of parent sources
close(h_parent); % Close the current waitbar box displaying progress through parent sources

for ps=1:n_sources % For each parent source
    data=H(:,:,ps); % Rename impulse response as "data" for saving as Roomsim audio MAT file. 
    save([H_filename '_S' num2str(ps)],'Fs','data'); % Save the sampling frequency and one L&R pair per data file..
end;% of ps counter loop for number of parent sources

%------------ Trim large arrays to remove trailing zeros ------------------
trim_length = max(len_source);
source = source(:,1:trim_length,:);
p_isource = p_isource(1:trim_length,:,:);
%--------------------------------------------------------------------------

%---------------------------------- End of roomsim_core.m --------------------------------------------
