function [df,tval,pval,cd,f] = doTTest(op,dv,iv,ivVal,paired,tail)

if nargin < 6; tail   = "both"; end
if nargin < 5; paired = false; end
if nargin < 4; error('Please assign a IV values'); end
if nargin < 3; error('Please assign a IV'); end
if nargin < 2; error('Please assign a DV'); end

try
    d1 = op.(dv)(op.(iv) == ivVal{1});
    d2 = op.(dv)(op.(iv) == ivVal{2});
catch
    error('Format error!')
end

if paired
    [~,pval,~,stat] = ttest(d1,d2,"Tail",tail);
else
    [~,pval,~,stat] = ttest2(d1,d2,"Tail",tail);
end
df   = stat.df;
tval = stat.tstat;
cd   = cohens_d(d1,d2,paired);

fmt = 't(%.0f) = %.02f, p = %.03f, d = %.02f \n';
f = sprintf(fmt, df, tval, pval, cd);
fprintf(f)
end

function cd = cohens_d(d1,d2,paired)

    m1 = mean(d1);
    m2 = mean(d2);
    m_diff = m1 - m2;

    if paired
        cd = m_diff / std(d1 - d2);
    else
        n1 = length(d1);
        n2 = length(d2);
    
        s1 = std(d1, 1);
        s2 = std(d2, 1);
        sp = sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1+n2-2) );
    
        cd = m_diff / sp;
    end
end

