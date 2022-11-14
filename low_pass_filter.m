function [data] = low_pass_filter(fig_loc);
% Low-Pass 4th order Butterworth filter Cut-off 0.91*Nyquist frequency used for bandlimiting/anti-aliasing
%i.e. Fs=44.1 kHz Fc= 20kHz, Fs=10 kHz Fc= 4.5Khz
%-------------------------------------------------------
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
%------------------------------------------------------------
% Functions called:

%*** Global Declarations ***
global curr_fig; %Keep track of current figure number

file_spec={'*.wav','Windows PCM (*.wav)';'*.mat','MAT-files (*.mat)';'*.*','All Files (*.*)'}; % *.wav first as more likely
banner='Choose the file for filtering (*.mat or *.wav)';
[Fs, X]=readf(file_spec,banner); % Get the data file

curr_fig=get(0,'CurrentFigure'); %Get handle of current figure
if isempty(curr_fig),
    figure(1); %Start first figure
else
    figure(curr_fig+1); %Start new figure
end
set(gcf,'position',fig_loc); % Put figure at this location
clf; % Make sure it's a clean sheet
colordef white; % Force plot backgrounds to white

subplot(2,1,1);
plot([1:size(X,1)],X(:,1),'b-',[1:size(X,1)],X(:,2),'r:'); %Plot the data for both channels
ylabel('Amplitude');
V=axis;
axis([0 V(2) V(3) V(4)]);
title('Original Data');

% Low-Pass filter
load LP_coeffs;% Numerator & Denominator coefficients generated using [B,A]=butter(4,20000/22050)
data(:,1) = filter(B,A,X(:,1)); % Filter left channel
data(:,2) = filter(B,A,X(:,2)); % Filter right channel

subplot(2,1,2);
plot([1:size(data,1)],data(:,1),'b-',[1:size(data,1)],data(:,2),'r:'); %Plot the data for both channels
ylabel('Amplitude');
axis([0 V(2) V(3) V(4)]);
title('Low-pass Filtered Data');

drawnow; %Force completion of all previous figure drawing before continuing
banner='Save the filtered result as *.wav or *.mat';
file_spec={'*.wav','Windows PCM (*.wav)';'*.mat','MAT-files (*.mat)';'*.*','All Files (*.*)'}; % wav first as more likely
savef(file_spec,Fs,data,banner); % Allow save of .mat or .wav file types

%---------- End of low_pass_filter.m -----------------------------