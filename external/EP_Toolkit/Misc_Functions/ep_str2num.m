function outVar=ep_str2num(inVar)
%  outVar=ep_str2num(inVar)
%       Implements original(?) Matlab behavior of str2num wherein it passes numeric variables without throwing an error.
%       Additionally adds restriction parameter to avoid evaluating command inputs if Matlab is version R2022a or newer.
%
%Inputs:
%  inVar         : Incoming variable or either numeric or character type.
%Outputs:
%  outVar        : Variable converted to numeric variable.

%History:
%  by Joseph Dien (6/14/25)
%  jdien07@mac.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

outVar=inVar;

MATLABver=ver('MATLAB');
[MatVer, MatRel]=strtok(MATLABver.Version,'.');
MatRel=MatRel(2:end);

if ischar(inVar)
    if ((str2double(MatVer)>22) || ((str2double(MatVer)==22) && (str2double(MatRel) >=12)))
        outVar=str2num(inVar,Evaluation="restricted");
    else
        outVar=str2num(inVar);
    end
end