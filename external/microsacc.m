function [sac, radius] = microsacc(x,vel,VFAC,MINDUR)
        %-------------------------------------------------------------------
        %  FUNCTION microsacc.m
        %  Detection of monocular candidates for microsaccades;
        %
        %  INPUT:
        %   x(:,1:2)         position vector
        %   vel(:,1:2)       velocity vector
        %   VFAC             relative velocity threshold
        %   MINDUR           minimal saccade duration
        %
        %  OUTPUT:
        %   radius         threshold velocity (x,y) used to distinguish microsaccs
        %   sac(1:num,1)   onset of saccade
        %   sac(1:num,2)   end of saccade
        %   sac(1:num,3)   peak velocity of saccade (vpeak)
        %   sac(1:num,4)   horizontal component     (dx)
        %   sac(1:num,5)   vertical component       (dy)
        %   sac(1:num,6)   horizontal amplitude     (dX)
        %   sac(1:num,7)   vertical amplitude       (dY)
        %---------------------------------------------------------------------
        % SDS... VFAC (relative velocity threshold) E&M 2006 use a value of VFAC=5

        % compute threshold
        % SDS... this is sqrt[median(x^2) - (median x)^2]
        msdx = sqrt( median(vel(:,1).^2,'omitnan') - (median(vel(:,1),'omitnan'))^2 );
        msdy = sqrt( median(vel(:,2).^2,'omitnan') - (median(vel(:,2),'omitnan'))^2 );
        if msdx<realmin
            msdx = sqrt( mean(vel(:,1).^2,'omitnan') - (mean(vel(:,1),'omitnan'))^2 );
            if msdx<realmin
                disp(['TRIAL: ' num2str(jj) ' msdx<realmin in eyelinkAnalysis.microsacc']);
            end
        end
        if msdy<realmin
            msdy = sqrt( mean(vel(:,2).^2,'omitnan') - (mean(vel(:,2),'omitnan'))^2 );
            if msdy<realmin
                disp(['TRIAL: ' num2str(jj) ' msdy<realmin in eyelinkAnalysis.microsacc']);
            end
        end
        radiusx = VFAC*msdx;
        radiusy = VFAC*msdy;
        radius = [radiusx radiusy];

        % compute test criterion: ellipse equation
        test = (vel(:,1)/radiusx).^2 + (vel(:,2)/radiusy).^2;
        indx = find(test>1);

        % determine saccades
        % SDS..  this loop reads through the index of above-threshold velocities,
        %        storing the beginning and end of each period (i.e. each saccade)
        %        as the position in the overall time series of data submitted
        %        to the analysis
        N = length(indx);
        sac = [];
        nsac = 0;
        dur = 1;
        a = 1;
        k = 1;
        while k<N
            if indx(k+1)-indx(k)==1     % looks 1 instant ahead of current instant
                dur = dur + 1;
            else
                if dur>=MINDUR
                    nsac = nsac + 1;
                    b = k;             % hence b is the last instant of the consecutive series constituting a microsaccade
                    sac(nsac,:) = [indx(a) indx(b)];
                end
                a = k+1;
                dur = 1;
            end
            k = k + 1;
        end

        % check for minimum duration
        % SDS.. this just deals with the final set of above threshold
        %       velocities; adds it to the list if the duration is long enough
        if dur>=MINDUR
            nsac = nsac + 1;
            b = k;
            sac(nsac,:) = [indx(a) indx(b)];
        end

        % compute peak velocity, horizonal and vertical components
        for s=1:nsac
            % onset and offset
            a = sac(s,1);
            b = sac(s,2);
            % saccade peak velocity (vpeak)
            vpeak = max( sqrt( vel(a:b,1).^2 + vel(a:b,2).^2 ) );
            sac(s,3) = vpeak;
            % saccade vector (dx,dy)            SDS..  this is the difference between initial and final positions
            dx = x(b,1)-x(a,1);
            dy = x(b,2)-x(a,2);
            sac(s,4) = dx;
            sac(s,5) = dy;

            % saccade amplitude (dX,dY)         SDS.. this is the difference between max and min positions over the excursion of the msac
            i = sac(s,1):sac(s,2);
            [minx, ix1] = min(x(i,1));              %       dX > 0 signifies rightward  (if ix2 > ix1)
            [maxx, ix2] = max(x(i,1));              %       dX < 0 signifies  leftward  (if ix2 < ix1)
            [miny, iy1] = min(x(i,2));              %       dY > 0 signifies    upward  (if iy2 > iy1)
            [maxy, iy2] = max(x(i,2));              %       dY < 0 signifies  downward  (if iy2 < iy1)
            dX = sign(ix2-ix1)*(maxx-minx);
            dY = sign(iy2-iy1)*(maxy-miny);
            sac(s,6:7) = [dX dY];
        end
end