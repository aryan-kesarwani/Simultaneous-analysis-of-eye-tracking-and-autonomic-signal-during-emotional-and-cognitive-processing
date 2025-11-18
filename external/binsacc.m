    function [sac, monol, monor] = binsacc(sacl,sacr)
        %-------------------------------------------------------------------
        %  FUNCTION binsacc.m
        %
        %  INPUT: saccade matrices from FUNCTION microsacc.m
        %   sacl(:,1:7)       microsaccades detected from left eye
        %   sacr(:,1:7)       microsaccades detected from right eye
        %
        %  OUTPUT:
        %   sac(:,1:14)       binocular microsaccades (right eye/left eye)
        %   monol(:,1:7)      monocular microsaccades of the left eye
        %   monor(:,1:7)      monocular microsaccades of the right eye
        %---------------------------------------------------------------------
        % SDS.. The aim of this routine is to pair up msaccs in L & R eyes that are
        %       coincident in time. Some msaccs in one eye may not have a matching
        %       msacc in the other; the code also seems to allow for a msacc in one
        %       eye matching 2 events in the other eye - in which case the
        %       larger amplitude one is selected, and the other discarded.

        if size(sacr,1)*size(sacl,1)>0

            % determine saccade clusters
            TR = max(sacr(:,2));
            TL = max(sacl(:,2));
            T = max([TL TR]);
            s = zeros(1,T+1);
            for i=1:size(sacl,1)
                s(sacl(i,1)+1:sacl(i,2)) = 1;   % SDS.. creates time-series with 1 for duration of left eye  msacc events and 0 for duration of intervals
                % NB.   sacl(i,1)+1    the +1 is necessary for the diff function [line 219] to correctly pick out the start instants of msaccs
            end
            for i=1:size(sacr,1)
                s(sacr(i,1)+1:sacr(i,2)) = 1;   % SDS.. superimposes similar for right eye; hence a time-series of binocular events
            end                                 %   ... 'binocular' means that either L, R or both eyes moved
            s(1) = 0;
            s(end) = 0;
            m = find(diff(s~=0));   % SDS.. finds time series positions of start and ends of (binocular) m.saccd phases
            N = length(m)/2;        % SDS.. N = number of microsaccades
            m = reshape(m,2,N)';    % SDS.. col 1 is all sacc onsets; col 2 is all sacc end points

            % determine binocular saccades
            NB = 0;
            NR = 0;
            NL = 0;
            sac = [];
            monol = [];
            monor = [];
            % SDS..  the loop counts through each position in the binoc list
            for i=1:N                                               % ..  'find' operates on the sacl (& sacr) matrices;
                l = find( m(i,1)<=sacl(:,1) & sacl(:,2)<=m(i,2) );  % ..   finds position of msacc in L eye list to match the timing of each msacc in binoc list (as represented by 'm')
                r = find( m(i,1)<=sacr(:,1) & sacr(:,2)<=m(i,2) );  % ..   finds position of msacc in R eye list ...
                % ..   N.B. some 'binoc' msaccs will not match one or other of the monoc lists
                if length(l)*length(r)>0                            % SDS..   Selects binoc msaccs.  [use of 'length' function is a bit quaint..  l and r should not be vectors, but single values..?]
                    ampr = sqrt(sacr(r,6).^2+sacr(r,7).^2);         % ..      is allowing for 2 (or more) discrete monocular msaccs coinciding with a single event in the 'binoc' list
                    ampl = sqrt(sacl(l,6).^2+sacl(l,7).^2);
                    [h ir] = max(ampr);                             % hence r(ir) in L241 is the position in sacr of the larger amplitude saccade (if there are 2 or more that occurence of binoc saccade)
                    [h il] = max(ampl);                             % hence l(il) in L241 is the position in sacl of the larger amplitude saccade (if there are 2 or more that occurence of binoc saccade)
                    NB = NB + 1;
                    sac(NB,:) = [sacr(r(ir),:) sacl(l(il),:)];      % ..      the final compilation selects the larger amplitude msacc to represent the msacc in that eye
                else
                    % determine monocular saccades
                    if isempty(l)                                 % If no msacc in L eye
                        NR = NR + 1;
                        monor(NR,:) = sacr(r,:);                    %..  record R eye monoc msacc.
                    end
                    if isempty(r)                                 %If no msacc in R eye
                        NL = NL + 1;
                        monol(NL,:) = sacl(l,:);                    %..  record L eye monoc msacc
                    end
                end
            end
        else
            % special cases of exclusively monocular saccades
            if size(sacr,1)==0
                sac = [];
                monor = [];
                monol = sacl;
            end
            if size(sacl,1)==0
                sac = [];
                monol = [];
                monor = sacr;
            end
        end
    end