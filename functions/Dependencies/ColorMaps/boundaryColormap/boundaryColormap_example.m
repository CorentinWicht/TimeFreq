clearvars
close all

x = linspace(-2, 2, 500);
y = linspace(-2, 3, 500);
[X, Y] = meshgrid(x, y);
Z = peaks(X, Y);

map = boundaryColormap(2, Z, 20);
imagesc([-2 2], [-2 3], Z)
set(gca, 'YDir', 'normal')
colormap(map)
colorbar