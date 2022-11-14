function roomsim_old2new
%Usage: roomsim_old2new;   Convert filename1.mat to filename2.wav or vice versa.
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

file_spec={'*.mat','MAT-files (*.mat)'}; %mat only
banner='Choose the OLD *.mat impulse response file to convert';

%------------------Get the file to convert--------------------
dot_ext=[];
while isempty(dot_ext),
    [name path] = uigetfile(file_spec, banner); %Display the dialogue box
    if ~any(name), 
        return;%**Alternate return for cancel operation. No data read from file.
    end;
    filename = lower([path name]); % Force path+filename to lowercase
    dot_ext = findstr(filename,'.'); % Find the extension marker
    if isempty(dot_ext) | length(filename)~=dot_ext+3 | ~strcmp(filename(dot_ext+1:dot_ext+3),'mat')
        h=errordlg('Only files with 3 letter *.mat extension allowed','roomsim_old2new error');
        beep;
        uiwait(h);
        dot_ext=[]; % Force resubmit of filename
    end;
end; % Of while loop to trap non .mat or .wav files
switch filename(dot_ext+1:dot_ext+3),
        
    case 'mat'
        Hp=[]; Fs=[]; % These are the named variables saved in a Roomsim audio *.mat data file
        load(filename); % Load an OLD ".mat" file and its sampling frequency
        
    otherwise, %shouldn't ever reach here but just in case
        h=errordlg('Data file extension not recognised, exiting','readf error');
        beep;
        uiwait(h);
        return;
        
end;
%--------------------------------------------------------------------

if isempty(Hp)|isempty(Fs),
    message = 'File is not Roomsim old audio format (does not contain the variable Hp or Fs)';
    h=errordlg(message,'roomsim_old2new error');
    beep;
    uiwait(h);
    return; % *** Alternate return
end;
data=Hp; % convert the variable name
banner='Save file as *.mat';
file_spec={'*.mat','MAT-files (*.mat)'}; % mat only
savef(file_spec,Fs,data,banner); % Allow save of NEW .mat file
%-------------------------- End of roomsim_old2new.m --------------------------