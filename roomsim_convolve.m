function roomsim_convolve(fig_loc);
%Usage: roomsim_convolve;   Convolves two stereo files (of type .mat or .wav),
% plots them, and the result, then save the result to a file, and displays spectrograms.
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
% Functions called: readf.m, savef.m

%*** Global Declarations ***
global curr_fig; %Keep track of current figure number

file_spec={'*.mat','MAT-files (*.mat)';'*.wav','Windows PCM (*.wav)';'*.*','All Files (*.*)'}; %mat first as more likely
banner='Choose the impulse response file (*.mat or *.wav)'; % Usually the shorter of the two files (for speed)
[Fs1, data]=readf(file_spec,banner); % Get the impulse response

if isempty(data),
    message = 'File_1 is not Roomsim audio format (does not contain the variable "data")';
    h=errordlg(message,'roomsim_convolve error');
    beep;
    uiwait(h);
    return; % *** Alternate return
end;

Hp=data; % Rename as the impulse response

if isempty(Fs1),
    banner = 'Sampling frequency not specified in file';
    prompt = {'Enter Sampling frequency (Hz): '};
    lines = 1;
    def = {'44100'}; %Default value
    answer = inputdlg(prompt,banner,lines,def);
    Fs1=str2num(answer{1});
end;

file_spec={'*.wav','Windows PCM (*.wav)';'*.mat','MAT-files (*.mat)';'*.*','All Files (*.*)'}; %wav first as more likely
banner='Choose the signal file (*.wav or *.mat)';
[Fs2, data]=readf(file_spec,banner); % Get the signal to be convolved with the impulse response

if isempty(data),
    message = 'File_2 is not Roomsim audio format (does not contain the variable "data")';
    h=errordlg(message,'roomsim_convolve error');
    beep;
    uiwait(h);
    return; % *** Alternate return
end;

if isempty(Fs2),
    banner = 'Sampling frequency not specified in file';
    prompt = {'Enter Sampling frequency (Hz): '};
    lines = 1;
    def = {'44100'}; %Default value
    answer = inputdlg(prompt,banner,lines,def);
    Fs2=str2num(answer{1});
end;

if abs(Fs1-Fs2) > eps,     
    banner='Roomsim_convolve';
    message=strvcat('Sampling frequencies are not equal','Data will be plotted against sample index number');
    h=msgbox(message,banner,'warn');  %Warn & beep.
    beep;
    uiwait(h);% Wait for user to acknowledge
    t_unit = false;
    Fs=[];
else
    banner='Roomsim_convolve';
    message='Use time for line Plots'' x-axis?';
    button=[];
    while isempty(button) % Loop until valid choice has been made (disable close window button)
        button = questdlg(message,banner,'Yes','No','Yes');
        beep;
    end; % Of while loop to disable inappropriate close window button operation
    if strcmp(button,'Yes') % 
        t_unit = true;
    else
        t_unit = false;
    end;
end;

%----- Do the convolutions using the FILTER primitive for speed. Hp is impulse response, data is (longer) audio file ----
data=[data; zeros(size(Hp,1)-1,2)]; % Append filter_length-1 zeros so convolution length is data_length+impulse_length-1
result=zeros(size(data));

msg_title = 'roomsim_convolve';
wait_msg_con ='Convolving Impulse Response and data please wait';
h_conv = waitbar(0.1,wait_msg_con,'name',msg_title); % Let the user know something is happening
result(:,1)=filter(Hp(:,1), 1, data(:,1)); % Convolve Left channel signals 
waitbar(0.5,h_conv); % Let the user know something is happening
result(:,2)=filter(Hp(:,2), 1, data(:,2)); % Convolve Right channel signals
waitbar(0.9,h_conv); % Let the user know something is happening
close(h_conv);

%-----------------------Plot the impulse response, data and convolution result for each channel---------------
%Set up the requested x axis scale,
conv_len=size(result,1); % Set to length of convolved file to align time plots
y_max_result = max(max(abs(result))); % Get the largest amplitude for vertical scaling of plots

if (t_unit == false)|(abs(Fs1-Fs2) > eps ), % sample number
    Fs=Fs2; % Use sampling frequency for data
    T=1/Fs; % Sampling period
    xscale1=[0:(size(Hp,1)-1)];
    xscale2=[0:(size(result,1)-1)]; 
    xtext=['Sample Index Number'];
else % time
    Fs=Fs1;
    T=1/Fs; % Sampling period
    xscale1=[0:(size(Hp,1)-1)]'*T;
    xscale2=[0:(size(result,1)-1)]'*T;
    xtext=['Time (s)'];
end  

curr_fig=get(0,'CurrentFigure'); %Get handle of current figure
if isempty(curr_fig),
    figure(1); %Start first figure
else
    figure(curr_fig+1); %Start new figure
end
set(gcf,'position',fig_loc); % Put figure at this location
clf ; %make sure it's a clean sheet
colordef white; % Force plot backgrounds to white

hold on
ytext=['Amplitude'];
Hp=[Hp; zeros(size(result,1)-size(Hp,1),2)]; % Extend the impulse response to same length as covolution result
subplot(3,1,1);
plot(xscale2,Hp(:,1),'b-',xscale2,Hp(:,2),'r:'); %Plot the impulse response of both channels L blue, R red
legend('Left Channel (blue)','Right Channel (red)');
V=axis;
axis([0 V(2) V(3) V(4)]);
title('Impulse response');
ylabel(ytext);

subplot(3,1,2);
plot(xscale2,data(:,1),'b-',xscale2,data(:,2),'r:'); %Plot the data for both channels
ylabel(ytext);
V=axis;
axis([0 V(2) -y_max_result y_max_result]);
title('Audio Data');

subplot(3,1,3);
plot(xscale2,result(:,1),'b-',xscale2,result(:,2),'r:'); %Plot the result for both channels
xlabel(xtext);
ylabel(ytext);
axis([0 V(2) -y_max_result y_max_result]);
title('Convolution result');

%---------------- Display spectrogram page --------------------------

% Normalise data and result for the spectrogram displays (removes relative constant gain factor)
[DD1,F1,Tf1] = specgram(data(:,1),[],Fs,[],[]); % Display data for spectrogram1
inv_max_DD1 = 1/max(max(abs(DD1)));
[DD2,F2,Tf2] = specgram(data(:,2),[],Fs,[],[]); % Display data for spectrogram2
inv_max_DD2 = 1/max(max(abs(DD2)));
[DD3,F3,Tf3] = specgram(result(:,1),[],Fs,[],[]); % Display data for spectrogram3
inv_max_DD3 = 1/max(max(abs(DD3)));
[DD4,F4,Tf4] = specgram(result(:,2),[],Fs,[],[]); % Display data for spectrogram4
inv_max_DD4 = 1/max(max(abs(DD4)));
inv_max_DD=min([inv_max_DD1 inv_max_DD2 inv_max_DD3 inv_max_DD4]); % Select an appropriate global scale factor
dB_DD1=20*log10(abs(DD1).*inv_max_DD+eps); % Convert to dB below 0dB
dB_DD2=20*log10(abs(DD2).*inv_max_DD+eps);
dB_DD3=20*log10(abs(DD3).*inv_max_DD+eps);
dB_DD4=20*log10(abs(DD4).*inv_max_DD+eps);

curr_fig=get(0,'CurrentFigure'); %Get handle of current figure
if isempty(curr_fig),
    figure(1); %Start first figure
else
    figure(curr_fig+1); %Start new figure
end
set(gcf,'position',fig_loc); % Put figure at bottom left of screen
clf ; %make sure it's a clean sheet
colormap jet;
clims=[-90 0]; % Set dB limits for spectrogram colormap
xtext=['Time (s)'];

subplot(2,2,1);
imagesc(Tf1,F1,dB_DD1,clims);
axis xy;
V=axis; % Get axis scaling for all spectrogram plots
axis([0 V(2) 0 V(4)]); 
ytext1=['Frequency (Hz)'];
title('Audio Data L');
xlabel([]); % No label
ylabel(ytext1);

subplot(2,2,2);
imagesc(Tf2,F2,dB_DD2,clims);
axis xy;
axis([0 V(2) 0 V(4)]);
title('Audio Data R');
xlabel([]); % No label
ylabel([]); % No label
colorbar;

subplot(2,2,3);
imagesc(Tf3,F3,dB_DD3,clims);
axis xy;
axis([0 V(2) 0 V(4)]);
title('Convolution result L');
xlabel(xtext);
ylabel(ytext1);

subplot(2,2,4);
imagesc(Tf4,F4,dB_DD4,clims);
axis xy;
axis([0 V(2) 0 V(4)]);
title('Convolution result R');
xlabel([xtext]);
ylabel([]); % No label
colorbar;
%--------- End plot spectrograms --------------------------
drawnow; %Force completion of all previous figure drawing before continuing

banner='To play the audio result, save as *.wav and use your media player';
file_spec={'*.wav','Windows PCM (*.wav)';'*.mat','MAT-files (*.mat)';'*.*','All Files (*.*)'}; % wav first as more likely
savef(file_spec,Fs,result,banner); % Allow save of .wav or .mat file types
%-------------------------- End of Roomsim_convolve.m --------------------------
