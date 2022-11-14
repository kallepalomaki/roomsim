function roomsim_convert
%Usage: roomsim_convert;   Convert filename1.mat to filename2.wav or vice versa.
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
%----------------------------------------------------------------------
% Functions called: readf.m, savef.m

file_spec={'*.mat','MAT-files (*.mat)';'*.wav','Windows PCM (*.wav)';'*.*','All Files (*.*)'}; %mat first as more likely
banner='Choose the file to convert (*.mat or *.wav)';
[Fs, data]=readf(file_spec,banner); % Get the file to convert

if isempty(data),
    message = 'File is not Roomsim audio format (does not contain the variable "data")';
    h=errordlg(message,'roomsim_convert error');
    beep;
    uiwait(h);
    return; % *** Alternate return
end;

if isempty(Fs),
    banner = 'Sampling frequency not specified in file';
    prompt = {'Enter Sampling frequency (Hz): '};
    lines = 1;
    def = {'44100'}; %Default value
    answer = inputdlg(prompt,banner,lines,def);
    Fs=str2num(answer{1});
end;

banner='Save file as *.wav or *.mat';
file_spec={'*.wav','Windows PCM (*.wav)';'*.mat','MAT-files (*.mat)';'*.*','All Files (*.*)'}; % wav first as more likely
savef(file_spec,Fs,data,banner); % Allow save of .wav or .mat file types
%-------------------------- End of roomsim_convert.m --------------------------