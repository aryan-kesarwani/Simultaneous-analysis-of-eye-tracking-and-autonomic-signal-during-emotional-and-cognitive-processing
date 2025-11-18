function op = getValues_sj(subjDir,dv,tw,iv,ivVal)
%%
% ----------------------------
% Get dependent variable matrix
% ----------------------------
dVar = getDV(subjDir, dv, tw);

% ----------------------------
% Multi-factor IV filtering
% ----------------------------
[labels,rmIdx] = getIV(subjDir,iv,ivVal);

% ----------------------------
% Compute mean value per condition
% ----------------------------
meanMat = nan(length(rmIdx), size(dVar,1));
num = nan(length(rmIdx),1);
for i = 1:length(rmIdx)
    meanMat(i,:) = mean(dVar(:,rmIdx{i}), 2, 'omitnan');
    num(i) = sum(rmIdx{i});
end

if ~isempty(tw)
    meanVal = meanMat;
else
    % meanVal = mean(meanMat, 2, 'omitnan');

    if strcmp(dv, 'Rate') == true
        meanVal = smoothdata(meanMat','movmean',100)';
    elseif ismember(dv, {'Amp','Vel','Ang','Cos','Sin'})
        meanVal = meanMat;
        if strcmp(dv,'Cos')
            meanMat = cos(meanMat);
        elseif strcmp(dv,'Sin')
            meanMat = sin(meanMat);
        end
        t = size(meanMat,2);
        for j = 1:size(meanMat,1)
            d_sub = meanMat(j,:);
            d_sub([1 end]) = mean(d_sub,'omitnan');
            t_sub = find(~isnan(d_sub));
            d_sub = d_sub(~isnan(d_sub));
            try
            d_sub = interp1(t_sub,d_sub,1:t,"linear");
            d_sub = smoothdata(d_sub,'movmean',500);
            meanVal(j,:) = d_sub;
            catch
            end
            
        end
    else
        meanVal = meanMat;
    end
end

% ----------------------------
% Output as table
% ----------------------------
op = [table(num, 'VariableNames', {'n'}), ...
    labels, table(meanVal, 'VariableNames', {dv})];