% This script creates an excel file containing the list of interpolated
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