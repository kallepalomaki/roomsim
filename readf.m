function [Fs, data] = readf(file_spec, title)
%Usage: [Fs, data] = readf(file_spec, title);   Read data from 
% filename.mat or filename.wav and place it in the array 'data'.
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
% Functions called: 

dot_ext=[];
while isempty(dot_ext),
    [name path] = uigetfile(file_spec, title); %Display the dialogue box
    if ~any(name), 
        return;%**Alternate return for cancel operation. No data read from file.
    end;
    filename = lower([path name]); % Force path+filename to lowercase
    dot_ext = findstr(filename,'.'); % Find the extension marker
    if isempty(dot_ext) | length(filename)~=dot_ext+3 | ...
            ( ~strcmp(filename(dot_ext+1:dot_ext+3),'wav') & ~strcmp(filename(dot_ext+1:dot_ext+3),'mat') )
        h=errordlg('Only files with 3 letter *.wav or *.mat extension allowed','readf error');
        beep;
        uiwait(h);
        dot_ext=[]; % Force resubmit of filename
    end;
end; % Of while loop to trap non .mat or .wav files
switch filename(dot_ext+1:dot_ext+3),
    case 'wav'
        [data, Fs, bits] = wavread(filename);	%Read a ".wav" file and load the sampling frequency
        
    case 'mat'
        data=[]; Fs=[]; % These are the named variables saved in a Roomsim audio *.mat data file
        load(filename); % Load a ".mat" file and keep the current sampling frequency
        
    otherwise, %shouldn't ever reach here but just in case
        h=errordlg('Data file extension not recognised, exiting','readf error');
        beep;
        uiwait(h);
        return;
        
end;

%---------------------------- End of readf.m ----------------------------------------