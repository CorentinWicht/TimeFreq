function Path = backslash (Loading_template)
%-------------------------------------------------------------------------%
% BACKSLASH
%-------------------------------------------------------------------------%
% This function doubles the number of slash '\' whenever it detects them in
% a string. This is sometimes useful as Matlab may crash for unknown
% reasons when setting up a path with only single slash.

% Usage:
%    >> Path = backslash (Loading_template);
%
% Inputs:
%   Loading_template       = String containing the path in which you want
%                            to double the '\'

% Outputs:
%   Path   - Returns the resulting string with doubled '\\'

% Author: Corentin Wicht, LCNS, 2018
% corentin.wicht@unifr.ch
    
j=1;
for i=1:length(Loading_template)
   if Loading_template(i)=='\'
       Path(j)=Loading_template(i);
       j=j+1;
       Path(j)='\';
   else
       Path(j)=Loading_template(i);
   end
   j=j+1;
end