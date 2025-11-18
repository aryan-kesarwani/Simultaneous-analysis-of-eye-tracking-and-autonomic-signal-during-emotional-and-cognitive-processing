function excludeOutliers(projectDir)
if nargin < 1
    projectDir = uigetdir(pwd, 'Please choose the project folder');
end


if ~ischar(projectDir) && ~isstring(projectDir)
    error('projectDir must be a string scalar or character vector.');
end

precDir = dir([projectDir '/derivatives/prec/']);


% Keep only directories with names matching 'sub-XX'
isSubFolder = startsWith({precDir.name}, 'sub-');
precDir = precDir(isSubFolder);

numSubj = length(precDir);

olFlag = true(numSubj,1);
for nSubj = 1:numSubj
    %%

    subIdx = sprintf('//sub-%02d',nSubj);
    subjDir = [precDir(nSubj).folder subIdx '//'];
    
    disp(['Subject Directory: ', subjDir]);

    subFile = dir([subjDir '/*beh.mat']);
    if isempty(subFile)
        warning('No behavior file found for subject %d', nSubj);
    end
    load([subFile.folder '//' subFile.name])
    
    disp(['Loaded behavior file: ', subFile.folder '//' subFile.name]);

    subFile = dir([subjDir '/*info.mat']);
    load([subFile.folder '//' subFile.name])

    subFile = dir([subjDir '/*rmFlag.mat']);
    load([subFile.folder '//' subFile.name])

    subFile = dir([subjDir '/*RateMatx.mat']);
    load([subFile.folder '//' subFile.name])

    ACC = logical(behData.ACC);
    rmFlag.ACC = ~ACC;

    rmFlag.RT = false(length(ACC),1);
    rmFlag.RT(ACC) = isoutlier(behData.RT(ACC),'mean');

    rmTable = [rmFlag.pup_l,rmFlag.ep,rmFlag.ACC,rmFlag.RT,rmFlag.sac];
    rmFlag.overall = any(rmTable')';

    subACC = [mean(behData.ACC(behData.cond==1)), ...
        mean(behData.ACC(behData.cond==2))];

    subInfo.olFlag = false;

    if any(subACC < 0.48)
        subInfo.olFlag = true;
    end

    if sum(RateMatx,"all") == 0
        subInfo.olFlag = true;
    end

    if subInfo.olFlag
        olFlag(nSubj) = false;
    end

    subjFileName = sprintf('sub-%02d_task-%s_',nSubj,subInfo.task);
    save([subjDir subjFileName 'info.mat'],       "subInfo")
    save([subjDir subjFileName 'rmFlag.mat'],     "rmFlag")
end
save([projectDir '/derivatives/prec/olFlag.mat'],       "olFlag")
