function Data = rmInsff(Data, num)

if size(table2array(Data(:,end)),2) > 1
    ol = isnan(mean(table2array(Data(:,end)),2));
else
    ol = isnan(table2array(Data(:,end)));
end
Data(ol,:) = [];

ol = unique(Data.Subj(Data.n < num));
Data(ismember(Data.Subj,ol),:) = [];

