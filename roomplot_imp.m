function roomplot_imp(Fs,H,sensor,fig_loc)
%Usage: roomplot_imp(Fs,H,sensor);
% Plot the impulse responses to each sensor
%------------------------------------------------------------------------------ 
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
%----------------------------------------------------------------------------------
% Functions called: 

%*** Global Declarations ***
global curr_fig; %Keep track of current figure number

%--------------------------- Initialisation ----------------------------
T=1/Fs; % Sampling period
n_plots=size(H,3); % Size finds number of pages ie sources (ps)

%Display the visualisation menu
menu_title='Roomplot_imp: Choose ordinate unit';
B1_text='Impulse response vs Time (seconds)';
B2_text='Impulse response vs Sample number';
beep;
M_IRLcase=0;
while M_IRLcase==0, % Loop until valid choice has been made (disable close window button)
    M_IRLcase = menu(menu_title,B1_text,B2_text);%Set up the requested x axis scale, sample number or time.
    switch M_IRLcase
        case 1 % Plot impulse response against time
            xscale=[0:size(H,1)-1]*T;
            xtext=['Time (s)'];
        case 2 % Plot impulse response against sample number
            xscale=[0:size(H,1)-1];
            xtext=['Sample Index Number'];
    end;
end; % Of while loop to disable inappropriate close window button operation

curr_fig=get(0,'CurrentFigure'); %Get handle of current figure
if isempty(curr_fig),
    figure(1); %Start first figure
else
    figure(curr_fig+1); %Start new figure
end
set(gcf,'position',fig_loc); % Put figure at this location
clf; % Make sure it's a clean sheet
colordef white; % Force plot backgrounds to white

for ps=1:n_plots %Plot the impulse responses for each original source
    subplot(n_plots,1,ps);
    hold on;
    if ps==1 
        title('Impulse response');
    end
    switch sensor % ---------- Identify receiver system---------------------
        case 'one_mic'
            plot(xscale,H(:,1,ps),'k-'); %Plot the impulse response to one sensor
            legend('Single Channel (black)');
        case {'two_mic','mithrir','cipicir'}  % Two sensor array or Head
            plot(xscale,H(:,1,ps),'b-',xscale,H(:,2,ps),'r:'); %Plot the impulse response to both channels
            switch sensor % ---------- Select receiver system---------------------
                case 'two_mic'
                    legend('Left Channel (blue)','Right Channel (red)');
                case {'mithrir','cipicir'}
                    legend('Left Ear (blue)','Right Ear (red)');
            end
        otherwise
            disp('Unknown receiver set up at impulse response plot');
            return
    end
    
    V=axis;
    axis([V(1) V(2) V(3) V(4)]);
    ytext=['Source ' num2str(ps)];
    xlabel(xtext); ylabel(ytext);
       
    hold off;
end;
drawnow; %Force completion of all previous figure drawing before continuing

%-------------------End of roomplot_imp.m ------------------------------------
