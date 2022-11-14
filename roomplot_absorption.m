function roomplot_absorption(F_abs,A,fig_loc)
% Usage: roomplot_absorption(F_abs,A,fig_loc);   Display the Frequency variation of the surface absorption coefficients
%-------------------------------------------------------------------------------------- 
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
%-----------------------------------------------------------------------------------------
% Functions called: 


%*** Global Declarations ***
global curr_fig; %Keep track of current figure number

curr_fig=get(0,'CurrentFigure'); %Get handle of current figure
if isempty(curr_fig),
    figure(1); %Start first figure
else
    figure(curr_fig+1); %Start new figure
end
set(gcf,'position',fig_loc); % Put figure at this locationn
clf; % Make sure it's a clean sheet
colordef white; % Force plot backgrounds to white
plot(F_abs,A(:,3),'r:x',F_abs,A(:,1),'m:o',F_abs,A(:,4),'b:*',F_abs,A(:,2),'g:s',F_abs,A(:,5),'k:^',F_abs,A(:,6),'y:v'); % Display the surface absorptions vs frequency
legend('Ay1','Ax1','Ay2','Ax2','Az1','Az2');
xlabel('Frequency (Hz)'); ylabel('Absorption Coefficient');
axis([0 max(F_abs) 0 1]);
title('Absorption Coefficients vs Frequency')
drawnow; %Force completion of all previous figure drawing before continuing

%-------------- End of roomplot_absorption.m ------------------------
