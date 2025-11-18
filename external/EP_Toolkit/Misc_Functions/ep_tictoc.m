function theDuration=ep_tictoc(theProc);
%  ep_tictoc(theProc);
%       Updates the status indicator.
%
%Inputs:
%  theProc      : String with current operation: 'load'=loading a file,'save'=saving a file,'done'=finished file operation.
%
%Outputs
%	theDuration: how long in seconds.

%History:
%  by Joseph Dien (6/11/19)
%  jdien07@mac.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 1999-2025  Joseph Dien
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global EPtictoc

theDuration=0;

if EPtictoc.stop==1
    theDuration=NaN;
    disp('Abort button pressed.')
    return;
end

if nargin==0
    theProc='tic';
end

try
    switch theProc
        case 'ioStart'
            set(EPtictoc.handles.reset,'BackgroundColor','blue');
            drawnow
        case 'ioFinish'
            set(EPtictoc.handles.reset,'BackgroundColor','cyan');
            drawnow
        case 'begin'
            set(EPtictoc.handles.reset,'BackgroundColor','cyan');
            drawnow
            EPtictoc.start=clock;
            EPtictoc.last=EPtictoc.start;
            EPtictoc.step=1;
        case 'end'
            set(EPtictoc.handles.reset,'BackgroundColor','green');
            set(EPtictoc.handles.T(EPtictoc.step),'BackgroundColor',[0.94 0.94 0.94]);
            drawnow
            EPtictoc.end=clock;
            theDuration=etime(EPtictoc.end,EPtictoc.start);
            if theDuration < 60
                fprintf('Process took %.2f seconds.\n',theDuration);
            else
                fprintf('Process took %0.2f minutes.\n',theDuration/60);
            end
        case 'tic'
            if isfield(EPtictoc,'last')
                if etime(clock,EPtictoc.last)>=1
                    EPtictoc.last=clock;
                    set(EPtictoc.handles.T(EPtictoc.step),'BackgroundColor',[0.94 0.94 0.94]);
                    EPtictoc.step=EPtictoc.step+1;
                    if EPtictoc.step==5
                        EPtictoc.step=1;
                    end
                    set(EPtictoc.handles.T(EPtictoc.step),'BackgroundColor','black');
                    drawnow
                end
            end
        case 'done'
            set(EPtictoc.handles.reset,'BackgroundColor',[0.94 0.94 0.94]);
            for iStep=1:4
                set(EPtictoc.handles.T(iStep),'BackgroundColor',[0.94 0.94 0.94]);
            end
            drawnow
        case 'error'
            set(EPtictoc.handles.reset,'BackgroundColor','red');
            drawnow
        otherwise
            disp('oops')
    end
catch ME
    msg{1}=['Aborting EPtictoc attempt to show progress on the process due to error:' ME.identifier];
    msg{2}=ME.message;
    if ~ishandle(EPtictoc.handles.reset)
        msg{3}='Note: Reset button missing.';
    end
    [msg]=ep_errorMsg(msg);
    EPtictoc.stop=1;
end



