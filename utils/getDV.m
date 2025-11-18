function dVar = getDV(subDir,dv,tw)

f = dir(fullfile(subDir,'*beh.mat'));
load(fullfile(f.folder,f.name))

f = dir(fullfile(subDir,'*info.mat'));
load(fullfile(f.folder,f.name))

if ismember(dv, behData.Properties.VariableNames)
    dVar = behData.(dv)';

elseif ismember(dv, fieldnames(subInfo))
    dVar = subInfo.(dv);

elseif strcmp(dv, 'Pupil_l')
    f = dir(fullfile(subDir,'*l_epbcData.mat'));
    load(fullfile(f.folder,f.name))
    dVar = l_epbcData;

elseif strcmp(dv, 'Pupil_r')
    f = dir(fullfile(subDir,'*r_epbcData.mat'));
    load(fullfile(f.folder,f.name))
    dVar = r_epbcData;

elseif strcmp(dv, 'tonicPupil_l')
    f = dir(fullfile(subDir,'*l_epData.mat'));
    load(fullfile(f.folder,f.name))
    dVar = l_epData;

elseif strcmp(dv, 'tonicPupil_r')
    f = dir(fullfile(subDir,'*r_epData.mat'));
    load(fullfile(f.folder,f.name))
    dVar = r_epData;

elseif strcmp(dv, 'Rate')
    f = dir(fullfile(subDir,'*RateMatx.mat'));
    load(fullfile(f.folder,f.name))
    dVar = RateMatx;
    if ~isempty(tw)
        dVar = dVar * (tw(2)-tw(1));
    end

elseif strcmp(dv, 'Amp')
    f = dir(fullfile(subDir,'*AmpMatx.mat'));
    load(fullfile(f.folder,f.name))
    dVar = AmpMatx;

elseif strcmp(dv, 'Vel')
    f = dir(fullfile(subDir,'*VelMatx.mat'));
    load(fullfile(f.folder,f.name))
    dVar = VelMatx;

elseif ismember(dv, {'Ang','Cos','Sin'})
    f = dir(fullfile(subDir,'*AngMatx.mat'));
    load(fullfile(f.folder,f.name))
    dVar = AngMatx;

else
    error('Unsupported DV: %s', dv);
end

if ~isempty(tw)
    tw = tw - subInfo.epoch(1);
    dVar = mean(dVar(tw(1):tw(2),:),'omitnan');
end