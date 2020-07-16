%% Author(s)

% Corentin Wicht (script, protocol) 
% GitHub : https://github.com/CorentinWicht

% If you have questions or want to contribute to this pipeline, feel free 
% to contact corentin.wicht@unifr.ch

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

%% This script creates an excel file containing the list of interpolated
% channels

% Empty cells
TempCell=cell(5,20);

% Copying data in the temporary cell matrix
for k=1:length(BadChannels_Parameters.Channels)
    for j=1:length(BadChannels_Parameters.Channels{k})
        TempCell(j,k)=BadChannels_Parameters.Channels{k}(j);
    end
end

% Creating a table
ExportTable=cell2table(TempCell);
ExportTable.Properties.VariableNames=cellfun(@(x) strrep(x,'\','_'),BadChannels_Parameters.Sbj,'UniformOutput',false);

% Write to excel file
writetable(ExportTable,'ListInterpChans.xlsx');