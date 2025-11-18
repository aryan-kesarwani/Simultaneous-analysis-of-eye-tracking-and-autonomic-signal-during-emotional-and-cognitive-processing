function [labels,rmIdx] = getIV(subjDir,iv,ivVal)
%%
% subFile = dir(fullfile(subDir,'*.mat'));
% for i = 1:length(subFile)
%     load(fullfile(subFile(i).folder,'\',subFile(i).name))
% end

f = dir(fullfile(subjDir,'*beh.mat'));
load(fullfile(f.folder,f.name))

f = dir(fullfile(subjDir,'*info.mat'));
load(fullfile(f.folder,f.name))

f = dir(fullfile(subjDir,'*rmFlag.mat'));
load(fullfile(f.folder,f.name))

numTrial = height(behData);

if isequal(iv, "none") || isempty(iv)
    rmIdx = mat2cell(~rmFlag.overall, numTrial, 1);
    labels = table();
else
    if isempty(ivVal)
        % Auto extract all unique combinations
        groupings = behData(:, iv);
        [~, ~, groupIdx] = unique(groupings, 'rows');
        nCond = max(groupIdx);
        rmIdx = cell(nCond,1);
        labels = unique(groupings);  % template row
        for i = 1:nCond
            mask = groupIdx == i;
            rmIdx{i} = ~rmFlag.overall & mask;
        end
    else
        % Use manually defined IV combinations
        rmIdx = cell(size(ivVal,1),1);
        for i = 1:size(ivVal,1)
            mask = true(numTrial,1);
            for j = 1:length(iv)
                mask = mask & behData.(iv{j}) == ivVal{i,j};
            end
            rmIdx{i} = ~rmFlag.overall & mask;
        end
        labels = cell2table(ivVal, 'VariableNames', iv);
    end
end



