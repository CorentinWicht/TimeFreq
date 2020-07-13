% NEED TO IMPROVE THE Level1, 2, folder comparison thing

function [Results,Cluster_Results]=Perm_Ttest(AllData,Design,EEG,varargin)

% This function computes (independent-, dependent-samples) t-tests relying 
% on permutation statistics + TFCE correction for multiple comparisons.
% It uses functions from the FMUT & ept_TFCE-MATLAB toolboxes.
% FMUT : https://github.com/ericcfields/FMUT
% ept_TFCE : https://github.com/Mensen/ept_TFCE-matlab

% Usage:
%    >> [Results,Pos,Cluster_Results]=Perm_Ttest(AllData,Design,EEG,...
%        varargin);
%
% Inputs (Mandatory):
%   ALLDATA = Structure containing each subject data in separate fields
%   (all conditions, groups, etc together)

%   DESIGN  = Design is a structure containing one field (Between or Within).
%   The field is a cell matrix of strings (Nx1) with the:
%   - 1st line containing either a 'B' (Between) or 'W' (Within)
%   - 2nd line containing the name of the factor (e.g. 'Group')
%   - 3rd to last lines containing the name of the two levels (e.g. OH, PBO)

%   EEG = EEGLab data structure (to retrieve the chanlocs and srate fields)

% Inputs (Optional):
%   N_PERMUTES  = [Integer] Number of permutations (def=5000)

%   PVAL = [Integer] Threshold for significance (def=0.05)

%   E_H = [Integer] TFCE correction parameters (def=[0.66 2])

%   ROOT FOLDER = [String] Upper folder path used to recompile the c-files 
%   using MEX, in case of errors. Will be bypassed if no argument is passed (def='').

% Outputs:
%   RESULTS = Data results structure containing the fields:
%   - F_Obs : Observed F-statistics
%   - TFCE_Obs : Observed TFCE-corrected F-statistics
%   - maxTFCE : Maximum TFCE values for each permutation maps 
%   - P_Values : Matrix of P-values

%   CLUSTER_RESULTS = Results of the clustering procedure 
%   (for details see function ept_calculateClusters)

% Author: Corentin Wicht, LCNS, 2019
% - corentin.wicht@unifr.ch
% - https://github.com/CorentinWicht


%% Set Defaults

N_Permutes  = 5000; % default number of permutations
Pval        = 0.05; % default threshold for significance
E_H         = [0.66 2]; % default parameters of E and H
root_folder = ''; % default upper path containing all files
 
% Process Secondary Arguments
if nargin > 4
    for i = 1:2:length(varargin)
        Param = varargin{i};
        Value = varargin{i+1};
        if ~ischar(Param)
          error('Flag arguments must be strings')
        end
        Param = lower(Param);
    
        switch Param
            case 'n_permutes'
                N_Permutes  = Value;
            case 'pval'
                Pval        = Value;
            case 'e_h'
                E_H         = Value;
            case 'root_folder'
                root_folder = Value;
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end

%% SETTINGS

% List of subjects data
SubjFN=fieldnames(AllData);

                        
% Determine if independent-/dependent-samples analysis
if isfield(Design, 'Between')
    CurrentDesign = Design.Between;
    StatsIdx = 'i';
elseif isfield(Design, 'Within')
    CurrentDesign = Design.Within;
    StatsIdx = 'd';
end

% Retrieving all the fields
AllFields = fieldnames(AllData);

% Groups/Conds to compare
CompareFN = CurrentDesign(3:end);

% For each level of the factor
for k=1:length(CompareFN)
    % Position of the current group/condition in the fields
    Idx = find(contains(AllFields,CurrentDesign(k+2)));
    % For each position
    for p=1:length(Idx)
        % Store the data in the big matrix
        Data{k}(:,:,p) = AllData.(AllFields{Idx(p)});
    end
    
    % Permuting data dimensions to match ept requirements
    Data{k} = permute(Data{k},[3 1 2]);
end

%% COMPUTING TEST STATISTICS
% For more information see:
% https://www.sciencedirect.com/science/article/pii/S1053811912010300?via%3Dihub

% Computing neighbouring channels 
ChN = ept_ChN2(EEG.chanlocs, 0);

% Errors might occur due to mex files
try
    % Running the permutation test
    % Info: flag_tf is for Time-Frequency (hence, the script will not look 
    % for channels (process them one by one)
    Results = ept_TFCE(Data{1},Data{2}, EEG.chanlocs,'nPerm', N_Permutes,...
        'rSample', EEG.srate,'ChN', ChN,'flag_tfce',1,'flag_tf',0,'type',...
        StatsIdx,'e_h',E_H);
catch
    if ~isempty(root_folder)
        % Errors may occur with mex files 
        % If error, recompile C-files manually

        % CD to the folder containing the ept_mex_TFCE.c files 
        TempPWD = pwd;
        cd(root_folder);
        EPTFolder = dir(['**/*' 'ept_TFCE-matlab-master/TFCE/Dependencies']);

        cd(EPTFolder.folder)
        mex -setup
        mex ept_mex_TFCE.c -compatibleArrayDims
        mex ept_mex_TFCE2D.c -compatibleArrayDims
        mex ept_mex_TFCE3D.c -compatibleArrayDims

        % CD back to the Data containing upper folder
        cd(TempPWD);

        % Running the permutation test
        Results = ept_TFCE(Data{1},Data{2}, EEG.chanlocs,'nPerm', N_Permutes,...
            'rSample', EEG.srate,'ChN', ChN,'flag_tfce',1,'flag_tf',0,'type',...
            StatsIdx,'e_h',E_H);
    else
        error('Consider providing a root_folder argument (path of upper folder) to recompile the c-files with MEX');
    end

end

%% Calculate clusters
try
    Cluster_Results = ept_calculateClusters(Results, ChN, Pval);
catch
    if ~isempty(root_folder)
        % CD to the folder containing the ept_ClustRes.c files 
        TempPWD = pwd;
        cd(root_folder);
        EPTFolder = dir(['**/*' 'ept_TFCE-matlab-master/ResultViewer/Dependencies']);

        cd(EPTFolder.folder)
        mex -setup
        mex ept_ClusRes.c -compatibleArrayDims
        mex ept_ClusRes3D.c -compatibleArrayDims
        
        Cluster_Results = ept_calculateClusters(Results, ChN, Pval);
        % CD back to the Data containing upper folder
        cd(TempPWD);
    else
        error('Consider providing a root_folder argument (path of upper folder) to recompile the c-files with MEX'); 
    end
end        

% if ispc 
%     % Plotting original and FDR-corrected results
%     % (useless)
%     Plot_test(EEG,CompareFN,CurrentDesign,Stats,Pval,Data,Normalization)
%

end