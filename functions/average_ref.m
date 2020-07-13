function EEG=average_ref(EEG,CzStruct,RefIndex)

% This function corrects the original pop_reref function flaws. Namely,
% pop_reref, replaces the old reference electrode at the end of the
% EEG.chanlocs structure. The current code replaces the old reference
% electrode at its right place (according to EEG.chanlocs.urchan)

% Usage:
%    >> [EEG] = average_ref(EEG,CzStruct);
%
% Inputs:
%   EEG       = EEG dataset structure
%   CzStruct  = EEG.chanlocs field of the old reference electrode
%   RefIndex  = Index of the old reference electrode from EEG.chanlocs.urchan
% Outputs:
%   EEG   - EEG dataset structure

% Author: Corentin Wicht, LCNS, 2018

    CzStruct.ref=CzStruct.labels;
    CzStruct.datachan=0;
    EEG = pop_reref( EEG, [],'refloc',CzStruct);
    [~,sortIndexes] = sort([EEG.chanlocs.urchan],'ascend');
    EEG.chanlocs = EEG.chanlocs(sortIndexes);
    EEG.data = EEG.data(sortIndexes,:,:);
end

