function retval=savef(file_spec, Fs, data, title)
%Usage: retval=savef(Fs, data, title);  Save the array 'data' to filename.mat or filename.wav.
%This function detects the file type extension and writes the desired format
%MATLAB loadable files are written if a .mat extension is detected
%WAV files are scaled into the range +/-1 then written as 16 bit with Fs 
%(The MATLAB warning message about clipping during wavwrite appears to be erroneous).
%------------------------------------------------------------------------------- 
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
%-----------------------------------------------------------------------------------
% Functions called: 

retval=1; [rows, cols] = size(data);	%Initialisation
if cols >2
    h=errordlg('savef can only handle mono or stereo files at present','savef error');
    beep;
    uiwait(h);
    return;
end

dot_ext=[];
while isempty(dot_ext),
    [name path] = uiputfile(file_spec, title); %Display the dialogue box
    if ~any(name), 
        return;%**Alternate return for cancel operation. No data saved to file.
    end;
    filename = lower([path name]); % Force path+filename to lowercase
    dot_ext = findstr(filename,'.'); % Find the extension marker
    if isempty(dot_ext) | length(filename)~=dot_ext+3 | ...
            ( ~strcmp(filename(dot_ext+1:dot_ext+3),'wav') & ~strcmp(filename(dot_ext+1:dot_ext+3),'mat') )
        h=errordlg('You must specify a 3 letter *.wav or *.mat extension','savef error');
        beep;
        uiwait(h);
        dot_ext=[]; % Force resubmit of filename
    end;
end; % Of while loop to trap missing 3 letter .extension

switch filename(dot_ext+1:dot_ext+3),
    case 'wav'
        max_data = max(max(abs(data))); %Find max of stereo data
        data=data/max_data;  %Scale into range +/- < 1
        
        if isempty(Fs),
            answer={};
            while isempty(answer),
                banner = 'Sampling frequency not specified for saving as *.wav';
                prompt = {'Enter a sampling frequency (Hz): '};
                lines = 1;
                def = {'8000'}; %Default value
                answer = inputdlg(prompt,banner,lines,def);
            end; % Of while loop to disable inappropriate CANCEL button operation
            Fs=str2num(answer{1});
        end;
        beep;
        
        wavwrite(data, Fs, 16, filename);	%Save as a ".wav" file, sampling frequency Fs, 16 bits
        
    case 'mat'
        % data is the name that Roomsim's readf.m will look for when reading a Roomsim audio *.mat file (avoids problem with compiled version)
        save(filename, 'Fs','data'); % Save as a Roomsim audio data ".mat" file
        
    otherwise, %shouldn't ever reach here but just in case
        h=errordlg('Data file extension not recognised, exiting','savef error');
        beep;
        uiwait(h);
        return;
        
end;
%--------------------------------- End of savef.m --------------------------------------
