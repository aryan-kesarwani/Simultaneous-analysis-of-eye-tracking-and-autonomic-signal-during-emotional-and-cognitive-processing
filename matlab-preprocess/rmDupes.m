function n_t = rmDupes(t)
% t = double(data{t_n}(1700:1725,1));
%%
% t = data(:,1);
t = double(t);
sep = find(diff(t)==0);
if ~isempty(sep)
    if sep(1) == 1
        t(sep(1)) = t(sep(1))-1;
        sep = find(diff(t)==0);
    end
end

for n = 1:length(sep)
    if (t(sep(n))-t(sep(n)-1)) >= 2
        t(sep(n)) = t(sep(n))-1;
    end
end

sep = find(diff(t)==0);
for n = 1:length(sep)-1
    if (t(sep(n)+2)-t(sep(n)+1)) == 2
        t(sep(n)+1) = t(sep(n)+1)+1;
    end
end

sep = find(diff(t)==0);
for n = 1:length(sep)
    t(sep(n)+1) = t(sep(n)+1)+1;
end
n_t = t;
