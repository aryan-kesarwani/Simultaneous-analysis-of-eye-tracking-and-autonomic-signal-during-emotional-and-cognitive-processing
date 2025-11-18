function [RateMatx, AmpMatx, VelMatx, AngMatx] = getSacMTX(rmSacData,tw)

numTrial  = length(rmSacData);
RateMatx  = zeros(tw(2)-tw(1),numTrial);
AmpMatx   = nan  (tw(2)-tw(1),numTrial);
VelMatx   = nan  (tw(2)-tw(1),numTrial);
AngMatx   = nan  (tw(2)-tw(1),numTrial);

for nTrial = 1:numTrial
    sac = rmSacData{nTrial};
    if isempty(sac)
        continue
    end
    sacIdx = sac(:,1)-tw(1)+1;
    sacAmp = sqrt(sac(:,6).^2 + sac(:,7).^2);
    sacVel = sac(:,3);
    sacAng = cart2pol(sac(:,6),sac(:,7));
    % sacAng(sacAng<0) = sacAng(sacAng<0) + 2*pi;

    RateMatx(sacIdx,nTrial) = 1;
    AmpMatx(sacIdx,nTrial) = sacAmp;
    VelMatx(sacIdx,nTrial) = sacVel;
    AngMatx(sacIdx,nTrial) = sacAng;
end
