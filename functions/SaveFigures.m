function SaveFigures(fig,Directory,Col,Formats)
%-------------------------------------------------------------------------%
% SAVING FULL SCREEN FIGURES
%-------------------------------------------------------------------------%
% This function enables to export full-screen pictures automically adapted
% to screen size. They will be exported in the directory specified by the
% argument Directory. Moreover, the functions can automatically save the
% figure in 3 different formats: pdf, bmp and eps. This function relies on
% another function called "backslash". 

% Usage:
%    >> [EEG] = SaveFigures(fig,Directory,Col,Formats);
%
% Inputs:
%   fig       = Figure to rescale and to save
%   Directory = Directory where to store the resulting figures 

% Optional input
%   Col       = Defines background color (e.g. 'b' for blue) (default is
%               white)
%   Formats   = Defines the format of export ('pdf', 'bmp' or 'eps')
%               (default is all three)


% Author: Corentin Wicht, LCNS, 2018
% corentin.wicht@unifr.ch

if nargin<3
  Col = 'w';
  Formats='All'; 
elseif nargin<4
  Formats='All';  
end

set(fig,'Units','Inches');
pos = get(fig,'Position');  
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(fig,'color',Col);
fig.InvertHardcopy = 'off';

    if strcmpi(Formats,'All') 
        % Saves the figure as as .pdf version
        print(fig,backslash(Directory),'-dpdf','-r0'); 
        % Saves the figure as .bmp format 
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(fig,backslash(Directory),'bmp');
        % Saves the topoplot as .eps format 
        saveas(fig,backslash(Directory),'epsc');
    elseif strcmpi(Formats,'pdf') 
        print(fig,backslash(Directory),'-dpdf','-r0'); 
    elseif strcmpi(Formats,'bmp')
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(fig,backslash(Directory),'bmp');
    elseif strcmpi(Formats,'eps') 
        % Saves the topoplot as .eps format 
        set(fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
        saveas(fig,backslash(Directory),'epsc');
    end

% Closes the current figure
close gcf
end
