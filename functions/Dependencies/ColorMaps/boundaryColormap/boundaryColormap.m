function map = boundaryColormap(boundary, data, N, map1, map2)
% Generates a colormap which would color the data below the given boundary
% using one colormap and color the data above the given boundary using
% another colormap.
% 
% The variable called 'boundary' is the target value above or belowe which
% a different colormap will be used;
% The variable called 'data' contains the data which has to be displayed
% using two colormaps. This variable can also just contain the maximum and
% minimum value of the data you would like to display with a colormap.
% The variable called 'N' gives a target to the size of the eventual
% colormap size (the number of colors in the colormap).
% The variables called 'map1' and 'map2' define functions, @(n) map1(n) and
% @(n) map2(n), which generates colormaps of length n.
% If 'map1' and 'map2' are not provided, then default colormaps will be
% used, which gives a reasonable contrast at the boundary edge and within
% the colors maps themself.

if ismember(nargin, [3 5])
    upper = max(data(:));
    lower = min(data(:));
    border = N / (upper - lower);
    size1 = round(border * (boundary - lower));
    size2 = round(border * (upper - boundary));
else
    assert(false, 'Incorrect number of inputs is used, either use 3 or 5 inputs.')
end

if nargin == 3
    map = [winter(size1)/2 + flipud(cool(size1))/2; flipud(autumn(size2))];
else
	map = [map1(size1); map2(size2)];
end

end