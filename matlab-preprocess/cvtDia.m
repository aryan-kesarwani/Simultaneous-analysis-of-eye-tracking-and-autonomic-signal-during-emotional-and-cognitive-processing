function nData = cvtDia(Data)

nData = Data;
p = polyfit([sqrt(6607) sqrt(3796) sqrt(1624) sqrt(414)],[8 6 4 2],1);

numTrial = length(Data);
for nTrial = 1:numTrial
    subData = Data{nTrial};
    subData(:,4) = double(polyval(p,sqrt(subData(:,4))));
    subData(:,9) = double(polyval(p,sqrt(subData(:,9))));
    nData{nTrial} = subData;
end


