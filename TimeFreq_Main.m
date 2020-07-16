%% Time-Frequency analyses on Epoched datasets

% NEED TO IMPLEMENT:
% 1) Correlations
% 2) Possibility to include more than 1 WS factors !
% 3) One-tailed tests ?
% 4) The export excel file should contain three columns/sheets when ANOVAs
% to account for the principal + interaction effects !

% BUGS :
% If folder comparison, might have a problem with SubjFN (i.e. same
% names!!)


% Version 0.4 / 13.07.2020 


%% ------------ XXX ---------- %%

%% Author(s)

% Corentin Wicht (script, protocol) 
% GitHub : https://github.com/CorentinWicht

% If you have questions or want to contribute to this pipeline, feel free 
% to contact corentin.wicht@unifr.ch

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

%% --------------------- PRESS F5 -------------------- %%
%% --------------------------------------------------- %%
clear variables
%% ----------------- PARAMETERS ------------- %%

% Get time
time_start = datestr(now);

% Get date + time
date_name = datestr(now,'dd-mm-yy_HHMM');

% ---------- SET DEPENDENCIES PATHS
addpath([pwd '/Functions']);
addpath([pwd '/Functions/Dependencies']);
addpath([pwd '/Functions/Dependencies/eeglab14_1_2b']);
addpath(genpath([pwd '/Functions/Dependencies/NMDv2-00']));
addpath(genpath([pwd '/Functions/Dependencies/ColorMaps']));
addpath(genpath([pwd '/Functions/Dependencies/ept_TFCE-matlab-master']));
addpath('Functions/Dependencies/FMUT-master');

% LOAD the parameters (from the prompt)
load([pwd '/TimeFreq_Param.mat'])

% ---------- PROMPTS
% Analyses 
evoked_ans = upper(answer{1});
induced_ans = upper(answer{2});
ITC_ans   = upper(answer{3});
avgref_ans = upper(answer{6});
Time_window = str2num(answer{7});
Normalization = upper(answer{8});

% Suffixes and parameters
SampRate = str2double(PromptInputs{1});
LowFreq = str2double(PromptInputs{2});  
HighFreq = str2double(PromptInputs{3});  
N_Permutes = str2double(PromptInputs{4});  
Pval = str2double(PromptInputs{5}); 
RelativeTolerance = str2double(PromptInputs{6});

% Frequency bands
BandsList=BandsList(all(~cellfun('isempty', BandsList),2),:);

% Path of the unique save folder (including date/time of analysis)
if ispc
    save_folder = [pwd '/Exports_' date_name];
    mkdir(['Exports_' date_name]);
    mkdir(['Exports_' date_name],'Raw');
    mkdir(['Exports_' date_name],'Stats');
else
    save_folder = [pwd '/Exports']; 
    mkdir('Exports');
    mkdir('Exports','Raw');
    mkdir('Exports','Stats');
end

% Loading the TF decomposed matrix
if strcmpi(ImportTF,'Yes')
    if ispc
        [file,path] = uigetfile([pwd '\Exports\Raw\*.mat'],'Select the file containing RAW data (RawTFData.mat in Raw folder)');
        Outputs = load([path file]);
    else
        Outputs = load([save_folder '/Raw/' 'RawTFData.mat']);
    end
end

% Path of most upper folder containing data
data_folder = uigetdir(pwd,'Select the most upper folder containing your .set EEG data files');
% data_folder = [pwd '/Data'];
root_folder = pwd;
cd(data_folder);

% List of files to analyze
FileList = dir(['**/*' '.set']);

% Some toolboxes included .set files in them
FileList(contains({FileList.name},'eeglab'))=[];

%% DATA REMOVAL (based on TimeFreq_Prompts.m)
% Filtering the data in FileList if a second factor was included
% Folders removal
if ~isempty(DesignList(:,3))
    DataRemove = DesignList(~cellfun('isempty',DesignList(:,3)),3);
    for m=1:length(DataRemove)
        FileList(contains({FileList.folder},DataRemove(m)))=[];
    end
end

% Files removal
if ~isempty(DesignList(:,4))
    DataRemove = DesignList(~cellfun('isempty',DesignList(:,4)),4);
    for m=1:length(DataRemove)
        FileList(contains({FileList.name},DataRemove(m)))=[];
    end
end

%% PARAMETERS SETTINGS
Parameters = [];
IdxParam = ones(1,length(FileList));
% Parameters
for j=1:length(FileList)
    AllParam = FileList(j).folder(length(data_folder)+2:end);
    if ispc
        AllParam = strsplit(AllParam,'\');
    else
        AllParam = strsplit(AllParam,'/');
    end

    % Looping through folders and storing them as levels in parameters
    for k=1:length(AllParam)
        Parameters.(sprintf('Level%d',k))(IdxParam(k),1)=AllParam(k);
        IdxParam(k) = IdxParam(k) + 1;
    end
end

% Fields of parameters (i.e. number of levels)
ParamFN = fieldnames(Parameters);

% Return error message if number of parameters > 3
if length(ParamFN) > 3
    error(['The number of data folders is to big (>3)! The current script can only accomodate 3 sub-folders (e.g. task, triggers, group).' ...
        newline 'Please reorder your files so that it matches the folder structure of the script.']) 
end

%% DESIGN

% Creating the Design structure separating Within and Between factors
Design=[];
for k=1:2 % i.e. 2 is the number of factors allowed
    if strcmpi(DesignList{1,k},'W')
        Design.Within=DesignList(~cellfun('isempty',DesignList(:,k)),k);
    elseif strcmpi(DesignList{1,k},'B')
        Design.Between=DesignList(~cellfun('isempty',DesignList(:,k)),k);
    end
end

% Create new exporting directories
mkdir([save_folder '/Raw'])
mkdir([save_folder '/Stats'])

%% ANALYSES

% Loading EEGLab functions
eeglab
close gcf

% set double-precision parameter
pop_editoptions('option_single', 0);

% Size of all the different files
sbj_high = length(FileList);

% Initialize variables
ProcessedCount = 0;

%% TF DECOMPOSITION

% Will execute if not imported (user defined)
if strcmpi(ImportTF,'No')
    
    % For each subject
    for sbj = 1:sbj_high

        % LOADING THE DATASET
        name_h = FileList(sbj).name(1:end-length('.set'));
        name_h(name_h == '_') = ' ';

        % Loading each file
        EEG = pop_loadset('filename',FileList(sbj).name,'filepath',FileList(sbj).folder);

        % Resampling (optional)
        if SampRate~=EEG.srate
            EEG = pop_resample(EEG, SampRate);
        end

        % Counting the number of files that can be processed (non-empty)
        if ~isempty(EEG.data)
           ProcessedCount = ProcessedCount + 1;
        end

        % INFOS
        fprintf('Currently computing time frequency decomposition for file %d out of %d\n',sbj,sbj_high);

        % Average referencing Cz
        if strcmp(avgref_ans,'Y')
            EEG = average_ref(EEG,EEG.chaninfo.nodatchans);
        end
        %% For each loaded/merged dataset
        
        ParamPath = FileList(j).folder(length(data_folder)+2:end);
        if ispc
            CurrentLevels = strsplit(ParamPath,'\')';
        else
            CurrentLevels = strsplit(ParamPath,'/')';
        end

        %% EVOKED ACTIVITY (after averaging over trials, phase-locked)

        if strcmpi(evoked_ans,'Y')

            % Averaging over the trials
            AvgEEG=EEG;
            AvgEEG.data=squeeze(mean(EEG.data,3));

            % First run to know the size of Freqs
            [~,ToDelete]=wt(AvgEEG.data(1,:),AvgEEG.srate,'fmin',LowFreq,'Padding',0,...
                'fmax',HighFreq,'Wavelet','Morlet','plot','off','Display','notify','f0',0.2,'RelTol',RelativeTolerance);

            % Reset the matrices
            TFdata=NaN([size(EEG.data,1),size(ToDelete,1),size(EEG.data,2)]);
            Freqs=NaN(size(ToDelete));
            clear ToDelete
            
            idx=1;
            textprogressbar('Calculating Evoked: ');

            % Computing Time-frequency (Wavelet Morlet) Decomposition
            % Channels X Freq X TimeFrames 
            % f0 = 0.2 is based on Schmiedt-Fehr et al., 2011a/2011b
            % Replaced Predictive by 0 padding, since predictive sometimes
            % freezed the analysis. 
            % For a comparison see: https://ars.els-cdn.com/content/image/1-s2.0-S1051200415000792-mmc1.pdf
            for l=1:size(EEG.data,1)
                [TFdata(l,:,:), Freqs]=wt(AvgEEG.data(l,:),AvgEEG.srate,'fmin',LowFreq,'Padding',0,...
                'fmax',HighFreq,'Wavelet','Morlet','plot','off','Display','notify','f0',0.2,'RelTol',RelativeTolerance);  
            
                % Display the progress
                textprogressbar((idx/(size(EEG.data,3)*size(EEG.data,1)))*100);
                idx=idx+1;
            end

            % Calculating Power of complex number (i.e. absolute value)
            % The structure depends on the number of parameters
            % Can go up to 3 levels/folders (e.g. task, triggers, group)
            if length(ParamFN) == 1
                % Saving the TF data in a structure
                Outputs.Evoked.(CurrentLevels{1}).Data.(strrep(name_h,' ','_'))=abs(TFdata);
                % Saving the frequency bins decomposition
                Outputs.Evoked.(CurrentLevels{1}).Freqs.(strrep(name_h,' ','_'))=Freqs;
            elseif length(ParamFN) == 2
                Outputs.Evoked.(CurrentLevels{1}).(CurrentLevels{2}).Data.(strrep(name_h,' ','_'))=abs(TFdata);
                Outputs.Evoked.(CurrentLevels{1}).(CurrentLevels{2}).Freqs.(strrep(name_h,' ','_'))=Freqs;
            elseif length(ParamFN) == 3
                Outputs.Evoked.(CurrentLevels{1}).(CurrentLevels{2}).(CurrentLevels{3}).Data.(strrep(name_h,' ','_'))=abs(TFdata);
                Outputs.Evoked.(CurrentLevels{1}).(CurrentLevels{2}).(CurrentLevels{3}).Freqs.(strrep(name_h,' ','_'))=Freqs;
            end
        end

        % Clearing the data matrix
        clear TFdata
        %% INDUCED ACTIVITY (before averaging over trials, non-phase locked)
        % aka EVENT-RELATED SPECTRAL POWER (ERSP)

        if strcmpi(induced_ans,'Y')

            % First run to know the size of Freqs
            [~,ToDelete]=wt(EEG.data(1,:,1),EEG.srate,'fmin',LowFreq,'Padding',0,...
                        'fmax',HighFreq,'Wavelet','Morlet','plot','off','Display','notify','f0',0.2);

            % Reset the matrices
            TFdata=NaN([size(EEG.data,1),size(EEG.data,3),size(ToDelete,1),size(EEG.data,2)]);
            Freqs=NaN(size(ToDelete));
            clear ToDelete
            idx=1;
            textprogressbar('Calculating Induced: ');

            % Looping over epochs
            for j=1:size(EEG.data,3)

                % Computing Time-frequency (Wavelet Morlet) Decomposition
                % Channels X Trials X Freq X TimeFrames
                for l=1:size(EEG.data,1)
                    % Calculating Power of complex number (i.e. absolute value)
                    [TFdata(l,j,:,:), Freqs]=wt(EEG.data(l,:,j),EEG.srate,'fmin',LowFreq,'Padding',0,...
                        'fmax',HighFreq,'Wavelet','Morlet','plot','off','Display','notify','f0',0.2);
                    
                   % Display the progress
                   textprogressbar((idx/(size(EEG.data,3)*size(EEG.data,1)))*100);
                   idx=idx+1;
                end
            end
            %% INTERTRIAL PHASE COHERENCE (ITC)
            % Based on fieldtrip
            % see : http://www.fieldtriptoolbox.org/faq/itc/
            if strcmpi(ITC_ans,'Y')
                               
                % Computing ITC
                ITC = TFdata./abs(TFdata); % divide by amplitude
                ITC = sum(ITC,2); % sum angles (trials)
                ITC = abs(ITC)/size(EEG.data,3); % take the absolute value and normalize
                ITC = squeeze(ITC); % remove the first singleton dimension
                                       
                % Saving the ITC data in a structure
                if length(ParamFN) == 1
                    % ITC magnitude is abs(itc); ITC phase in radians is angle(itc), or in deg phase(itc)*180/pi.
                    Outputs.ITC.(CurrentLevels{1}).Data.(strrep(name_h,' ','_'))=ITC;
                    % Saving the frequency bins decompositions
                    Outputs.ITC.(CurrentLevels{1}).Freqs.(strrep(name_h,' ','_'))=Freqs;
                elseif length(ParamFN) == 2
                    Outputs.ITC.(CurrentLevels{1}).(CurrentLevels{2}).Data.(strrep(name_h,' ','_'))=ITC;
                    Outputs.ITC.(CurrentLevels{1}).(CurrentLevels{2}).Freqs.(strrep(name_h,' ','_'))=Freqs;
                elseif length(ParamFN) == 3
                    Outputs.ITC.(CurrentLevels{1}).(CurrentLevels{2}).(CurrentLevels{3}).Data.(strrep(name_h,' ','_'))=ITC;
                    Outputs.ITC.(CurrentLevels{1}).(CurrentLevels{2}).(CurrentLevels{3}).Freqs.(strrep(name_h,' ','_'))=Freqs;
                end
            end

            % Averaging over the trials
            TFdata=squeeze(mean(abs(TFdata),2));
            
            if length(ParamFN) == 1
                % Saving the TF data in a structure
                Outputs.Induced.(CurrentLevels{1}).Data.(strrep(name_h,' ','_'))=TFdata;
                % Saving the frequency bins decomposition
                Outputs.Induced.(CurrentLevels{1}).Freqs.(strrep(name_h,' ','_'))=Freqs;
            elseif length(ParamFN) == 2
                Outputs.Induced.(CurrentLevels{1}).(CurrentLevels{2}).Data.(strrep(name_h,' ','_'))=TFdata;
                Outputs.Induced.(CurrentLevels{1}).(CurrentLevels{2}).Freqs.(strrep(name_h,' ','_'))=Freqs;
            elseif length(ParamFN) == 3
                Outputs.Induced.(CurrentLevels{1}).(CurrentLevels{2}).(CurrentLevels{3}).Data.(strrep(name_h,' ','_'))=TFdata;
                Outputs.Induced.(CurrentLevels{1}).(CurrentLevels{2}).(CurrentLevels{3}).Freqs.(strrep(name_h,' ','_'))=Freqs;
            end
        end
    end 
    
    % Deleting parameters duplicate
    ParamFN = fieldnames(Parameters);
    for k=1:numel(ParamFN)
        Parameters.(ParamFN{k})=unique(Parameters.(ParamFN{k}));
    end

    % Clearing unused variables
    if exist('TFdata','var')
        clear TFdata
    end
    if exist ('ITC','var') 
        clear ITC
    end

    % Exporting original values
    save([save_folder '/Raw/' 'RawTFData.mat'],'-struct','Outputs')
    
else
    % LOADING THE FIRST DATASET TO RETRIEVE THE LABELS
    name_h = FileList(1).name(1:end-length('.set'));
    name_h(name_h == '_') = ' ';

    % Loading each file
    EEG = pop_loadset('filename',FileList(1).name,'filepath',FileList(1).folder);

    % Resampling (optional)
    if SampRate~=EEG.srate
        EEG = pop_resample(EEG, SampRate);
    end

    % Average referencing Cz
    if strcmp(avgref_ans,'Y')
        EEG = average_ref(EEG,EEG.chaninfo.nodatchans);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------------------------------------%

% Decision about the statistical design 
DesignFN = fieldnames(Design);

% T-tests
if length(DesignFN) == 1
    if contains(DesignFN,'Within')
        StatTest = 'Dependent-samples t-test';
        StatsIdx = 'd';
    elseif contains(DesignFN,'Between')
        StatTest = 'Independent-samples t-test';
        StatsIdx = 'i';
    end
    
% ANOVA's
elseif length(DesignFN) == 2
    if sum(contains(DesignFN,'Within') + contains(DesignFN,'Between')) == 2
        StatTest = 'Mixed ANOVA';
        % StatsIdx = 'i';
    elseif sum(contains(DesignFN,'Within')) == 2
        StatTest = 'Repeated-Measures ANOVA';
    end
end

% FIELDNAMES %%%%%%%%%%%%%
% TF methods
AnalysesFN=fieldnames(Outputs);
% 1st Level (e.g. Task)
if length(ParamFN) == 1 || length(ParamFN) == 2 || length(ParamFN) == 3
    Level1FN=fieldnames(Outputs.(AnalysesFN{1}));
end
% 2nd Level (e.g. Triggers)
if length(ParamFN) == 2 || length(ParamFN) == 3
    Level2FN=fieldnames(Outputs.(AnalysesFN{1}).(Level1FN{1}));
else
    Level2FN = {''};
end
% % 3rd Level (e.g. Groups)
if length(ParamFN) == 3
    Level3FN=fieldnames(Outputs.(AnalysesFN{1}).(Level1FN{1}).(Level2FN{1}));
else
    Level3FN = {''};
end
% Subjects
SubjFN=cellfun(@(x) x(1:end-4),{FileList.name}','UniformOutput',0);

% If SUBJFN in different folders but with SAME name, will generate errors!!
% NEED TO FIND A WAY WHEN FOLDER COMPARISON! 

% Number of files processed (for the log)
if strcmpi(ImportTF,'Yes')
    ProcessedCount = length(SubjFN); 
end

% Preallocation
Loop = 1;
Num = 1;
DispText = '';
if contains(StatTest,'t-test')
    cluster_thresh=zeros(length(AnalysesFN)*length(Level1FN)*length(Level2FN)*length(Level3FN)*size(BandsList,1),1);
    NumSigData=zeros(length(AnalysesFN)*length(Level1FN)*length(Level2FN)*length(Level3FN)*size(BandsList,1),1);
elseif contains(StatTest,'ANOVA')
    cluster_thresh=zeros(length(AnalysesFN)*length(Level1FN)*length(Level2FN)*length(Level3FN)*size(BandsList,1),3);
    NumSigData=zeros(length(AnalysesFN)*length(Level1FN)*length(Level2FN)*length(Level3FN)*size(BandsList,1),3);
end

% For each Time-Frequency metrics
for k=1:length(AnalysesFN)
   
    % For each frequency band of interest
    for o=1:size(BandsList,1)
        
        % Looking for index of frequency boundaries
        CurrentBounds=str2num(BandsList{o,2});
    
        % For each 1st level
        for m=1:length(Level1FN)

            % For 2nd level
            for n=1:length(Level2FN)
                
                % For 3rd level
                for p=1:length(Level3FN)
                    
                    % Retrieving frequency band and data
                    if length(ParamFN) == 1
                        AllFreqs=Outputs.(AnalysesFN{k}).(Level1FN{m}).Freqs.(SubjFN{1});
                        AllData=Outputs.(AnalysesFN{k}).(Level1FN{m}).Data;
                        % Filling the Stats structure with file information
                        Stats.FileName=[AnalysesFN{k} '_' Level1FN{m} '_' BandsList{o,1}];
                    elseif length(ParamFN) == 2
                        AllFreqs=Outputs.(AnalysesFN{k}).(Level1FN{m}).(Level2FN{n}).Freqs.(SubjFN{1});
                        AllData=Outputs.(AnalysesFN{k}).(Level1FN{m}).(Level2FN{n}).Data;
                        Stats.FileName=[AnalysesFN{k} '_' Level1FN{m} '_' Level2FN{n} '_' BandsList{o,1}];
                    elseif length(ParamFN) == 3
                        AllFreqs=Outputs.(AnalysesFN{k}).(Level1FN{m}).(Level2FN{n}).(Level3FN{p}).Freqs.(SubjFN{1});
                        AllData=Outputs.(AnalysesFN{k}).(Level1FN{m}).(Level2FN{n}).(Level3FN{p}).Data;
                        Stats.FileName=[AnalysesFN{k} '_' Level1FN{m} '_' Level2FN{n} '_' Level3FN{p} '_' BandsList{o,1}];
                    end
                    
                    % Finding position of frequency boundaries index
                    BoundsIdx=find(ge(AllFreqs,CurrentBounds(1)) & le(AllFreqs,CurrentBounds(2)));
                    
                    % Only selecting data in frequency band of interest
                    AllData=structfun(@(x) squeeze(mean(x(:,BoundsIdx,:),2)),AllData,'UniformOutput',0);

                    % Restricting data length (optional)
                    if ~isempty(Time_window)
                        SavedEEGxmin=EEG.xmin;
                        IdxBeg = round(Time_window(1)*(EEG.srate/1000)) + abs(round(SavedEEGxmin*EEG.srate)); % second term is the baseline 
                        IdxEnd = round(Time_window(end)*(EEG.srate/1000)) + abs(round(SavedEEGxmin*EEG.srate));
                        AllData=structfun(@(x) x(:,IdxBeg:IdxEnd),AllData,'UniformOutput',0);

                        % Adjusting EEG structure
                        EEG.pnts = length(IdxBeg:IdxEnd);
                        EEG.times = EEG.times(IdxBeg:IdxEnd);
                        EEG.xmin = Time_window(1)/1000;
                        EEG.xmax = Time_window(end)/1000;
                    end

                    % Power spectrum normalization (optional)
                    if strcmpi(Normalization,'Y')
                        AllData=structfun(@(x) 10*log10(x),AllData,'UniformOutput',0);
                    end 
                    
                    % Calling the statistics scripts
                    if contains(StatTest,'t-test')
                        
                        % T-TEST
                        [Results,Cluster_Results]=Perm_Ttest(AllData,Design,...
                            EEG,'root_folder',root_folder,'N_Permutes',N_Permutes,...
                            'Pval',Pval);   
                        
                        % Permutation threshold (e.g. 95% confidence interval)                     
                        U = round((1-Pval)*N_Permutes); 
                        MaxTFCE=sort(Results.maxTFCE);
                        cluster_thresh(Loop)= MaxTFCE(U);
                        
                        % Determining number of significant results
                        NumSigData(Loop) = length(Cluster_Results);  
                        
                    elseif contains(StatTest,'ANOVA')
                        
                        % F-TEST
                        [Results,Cluster_Results]=Perm_ANOVA(AllData,Design,...
                            EEG,'root_folder',root_folder,'N_Permutes',N_Permutes,...
                            'Pval',Pval);
                        
                        % Permutation threshold (e.g. 95% confidence interval)
                        U = round((1-Pval)*N_Permutes); 
                        MaxTFCE.A=sort(Results.A.maxTFCE);
                        MaxTFCE.B=sort(Results.B.maxTFCE);
                        MaxTFCE.AB=sort(Results.AB.maxTFCE);
                        
                        % Storing the clustering threshold
                        cluster_thresh(Loop,1)= MaxTFCE.A(U);
                        cluster_thresh(Loop,2)= MaxTFCE.B(U);
                        cluster_thresh(Loop,3)= MaxTFCE.AB(U);
                        
                        % Determining number of significant results
                        NumSigData(Loop,1) = length(Cluster_Results.A); 
                        NumSigData(Loop,2) = length(Cluster_Results.B); 
                        NumSigData(Loop,3) = length(Cluster_Results.AB); 
                    end
                    
                    % FOLDER COMPARISON
                    % NEED TO TEST IF THIS IS ENOUGH !!! 
                    
                    % Ending the loop corresponding to the folder
                    % Otherwise the analysis will be run twice!
                    AllLevels = unique(cellfun(@(x) x(length(data_folder)+2:end),{FileList.folder},'UniformOutput',0))';
                    AllLevelsSplit = cellfun(@(x) strsplit(x,'\'),AllLevels, 'UniformOutput',0);
                    for u=1:length(AllLevelsSplit{:})
                        
                        % For each WS/BS factor(s)
                        DesignFields = fieldnames(Design);
                        for f=1:numel(DesignFields)
                            
                            % THIS WILL NOT WORK IF FACTORS DO NOT HAVE THE
                            % SAME LENGTH (i.e. number of levels) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                            for t=1:size(Design.(DesignFields{f}),2)
                                CurrentDesign = Design.(DesignFields{t});    
                                % Pos will indicate which level is
                                % concerned by the folder comparision
                                Pos = find(ismember(AllLevelsSplit{:}(u),CurrentDesign{end-1}));
                                if ~isempty(Pos)                 
                                    % Ending the corresponding loop
                                    % Will avoid to run twice the same analysis
                                    if Pos == 1
                                        m=length(Level1FN);
                                    elseif Pos == 2
                                        n=length(Level2FN);
                                    elseif Pos == 3
                                        p=length(Level3FN);       
                                    end 

                                    % Updating the stats filename (in case of folders
                                    % comparison)
                                    Stats.FileName=strrep(Stats.FileName,['_' CurrentDesign{end-1}],'');
                                    break
                                end
                            end
                        end
                    end

                    % Increase the loop counter
                    Loop=Loop+1;   
                    
                    % Storing results
                    Stats.TFCE = Results;  
                    Stats.Clusters = Cluster_Results;  
                    Stats.Test = StatTest;

                    % Exporting thresholded statistical values
                    save([save_folder '/Stats/' Stats.FileName '.mat'],'Stats')

                    % Free up some space for further analyses
                    clear Stats AllData TempData 
                end
            end
        end
    end
end

%% EXPORTS
SizeExp=length(AnalysesFN)*length(Level1FN)*length(Level2FN)*length(Level3FN)*size(BandsList,1);

% Exporting the results in an excel file
if contains(StatTest,'t-test')
    Export=cell(SizeExp, length(ParamFN)+3);
    if length(ParamFN) == 1
        Headers={'Analysis','FrequencyBand','Level1','ClusterThreshold','SigClusters'};
    elseif length(ParamFN) == 2
        Headers={'Analysis','FrequencyBand','Level1','Level2','ClusterThreshold','SigClusters'};
    elseif length(ParamFN) == 3
        Headers={'Analysis','FrequencyBand','Level1','Level2','Level3','ClusterThreshold','SigClusters'};
    end
elseif contains(StatTest,'ANOVA')
    Export=cell(SizeExp, length(ParamFN)+5); % For principal effects + interaction    
    if length(ParamFN) == 1
        Headers={'Analysis','FrequencyBand','Level1','ClusterThreshold_A','ClusterThreshold_B','ClusterThreshold_AB',...
            'SigClusters_A','SigClusters_B','SigClusters_AB'};
    elseif length(ParamFN) == 2
        Headers={'Analysis','FrequencyBand','Level1','Level2','ClusterThreshold_A','ClusterThreshold_B','ClusterThreshold_AB',...
            'SigClusters_A','SigClusters_B','SigClusters_AB'};
    elseif length(ParamFN) == 3
        Headers={'Analysis','FrequencyBand','Level1','Level2','Level3','ClusterThreshold_A','ClusterThreshold_B','ClusterThreshold_AB',...
            'SigClusters_A','SigClusters_B','SigClusters_AB'};
    end
end

% Regulate position inside the excel sheet (i.e. columns)

Pos = 1;
% Analysis
TemporDat=cellfun(@(x) repmat({x},[SizeExp/length(AnalysesFN),1]),AnalysesFN,'UniformOutput',false);
Export(:,Pos)=vertcat(TemporDat{:});
Pos = Pos + 1;

% FrequencyBand
% TemporDat=cellfun(@(x) repmat({x},[(length(AnalysesFN)),1]),...
%     BandsList(:,1),'UniformOutput',false);
TemporDat=repmat(BandsList(:,1),[length(AnalysesFN)*length(Level1FN)*length(Level2FN)*length(Level3FN),1]);
%Export(:,Pos)=vertcat(TemporDat{:});
Export(:,Pos)=TemporDat;
Pos = Pos + 1;

% Level 1
TemporDat=cellfun(@(x) repmat({x},[SizeExp/(length(Level1FN)*length(AnalysesFN)),1])...
    ,Level1FN,'UniformOutput',false);
TemporDat=repmat(TemporDat,[length(AnalysesFN),1]);
Export(:,Pos)=vertcat(TemporDat{:});
Pos = Pos + 1;

% Level 2
if length(ParamFN) > 1
    TemporDat=cellfun(@(x) repmat({x},[SizeExp/(length(Level2FN)*length(Level1FN)*length(AnalysesFN)),1]),...
        Level2FN,'UniformOutput',false);
    TemporDat=repmat(TemporDat,[length(AnalysesFN)*length(Level1FN),1]);
    Export(:,Pos)=vertcat(TemporDat{:});
    Pos = Pos + 1;
end

% Level 3
if length(ParamFN) > 2
    TemporDat=cellfun(@(x) repmat({x},[SizeExp/(length(Level3FN)*length(Level2FN)...
        *length(Level1FN)*length(AnalysesFN)),1]),Level3FN,'UniformOutput',false);
    TemporDat=repmat(TemporDat,[length(AnalysesFN)*length(Level1FN),1]);
    Export(:,Pos)=vertcat(TemporDat{:});
    Pos = Pos + 1;
end

% Cluster Thresholds
if contains(StatTest,'t-test')
    TemporDat=reshape(cluster_thresh',[size(cluster_thresh,1)*size(cluster_thresh,2),1]);
    Export(:,Pos)=num2cell(TemporDat);
    Pos = Pos + 1;  
elseif contains(StatTest,'ANOVA')
    Export(:,Pos:Pos+2)=num2cell(cluster_thresh);
    Pos = Pos + 3;    
end

% Number of Significant clusters
if contains(StatTest,'t-test')
    TemporDat=reshape(NumSigData',[size(NumSigData,1)*size(NumSigData,2),1]);
    Export(:,Pos)=num2cell(TemporDat);
elseif contains(StatTest,'ANOVA')
    Export(:,Pos:Pos+2)=num2cell(NumSigData); 
end

% Create the table
TableExp=cell2table(Export);
TableExp.Properties.VariableNames=Headers;

% Write to excel file
format long g
CurrentDateTime=fix(clock);
writetable(TableExp,[save_folder '/' sprintf('Analyses_%s_%dH%d.xlsx',datetime('today'),...
    CurrentDateTime(end-2),CurrentDateTime(end-1))]);

%% Log
time_end = datestr(now);
username=getenv('USERNAME');
Computer=computer;

% Creating the log file
date_name = datestr(now,'dd-mm-yy_HHMM');
fid = fopen([save_folder '/TimeFreqlog_' date_name '.txt'],'w');

% date, starting time, finished time, number of analyzed files
fprintf(fid,'%s\t%s\r\n',['Start : ',time_start],['End: ',time_end]);
fprintf(fid,'\r\n%s\r\n',['Username : ' username]);
fprintf(fid,'%s\r\n',['Computer type : ' Computer]);

if ~exist('TFMatrix','var')
    fprintf(fid,'\r\n%s',[num2str(ProcessedCount) ' file(s) were PROCESSED out of ' num2str(sbj_high) '.']); 
else
    fprintf(fid,'\r\n%s',[num2str(ProcessedCount) ' file(s) were IMPORTED out of ' num2str(sbj_high) '.']); 
end

% Summary of analysis parameters
fprintf(fid,'\r\n\r\n%s\r\n','------ ANALYSIS PARAMETERS SUMMARY ------');
fprintf(fid,'\r\n%s',sprintf('Sampling rate : %d Hz',SampRate));
fprintf(fid,'\r\n%s',sprintf('Frequency range : %d-%d Hz',LowFreq,HighFreq));
fprintf(fid,'\r\n%s',sprintf('Number of permutations : %d',N_Permutes));
fprintf(fid,'\r\n%s',sprintf('Threshold of significance (p-value) : %.3f',Pval));
fprintf(fid,'\r\n%s',sprintf('Relative tolerance for the cone of influence (default = 0.01) : %.3f',RelativeTolerance));
if strcmpi(avgref_ans,'Y')
    fprintf(fid,'\r\n%s',sprintf('%d Files were average referenced.',ProcessedCount)); 
else
    fprintf(fid,'\r\n%s',sprintf('%d Files were not average referenced.',ProcessedCount)); 
end

% List of Files for which significant clusters were identified
fprintf(fid,'\r\n\r\n%s\r\n','------ STATISTICAL SIGNIFICANCE SUMMARY ------');
fprintf(fid,'%s\r\n',sprintf('Test statistic : %s',StatTest));
if isfield(Design,'Between')
    fprintf(fid,'%s\r\n',sprintf('Between-subject factor : %s Vs %s',Design.Between{end-1},Design.Between{end}));
end
if isfield(Design,'Within')
    fprintf(fid,'%s\r\n',sprintf('Within-subject factor : %s Vs %s',Design.Within{end-1},Design.Within{end}));
end
fprintf(fid,'\r\n%s\r\n','You will find below the list of frequency bands and conditions for which significant differences were identified:');

for k=1:SizeExp
    for m=1:size(NumSigData,2)
        if NumSigData(k,m)~=0  
            if length(ParamFN) == 1
                fprintf(fid,'\r\n%s',[Export{k,1} '_' Export{k,2} '_' Export{k,3} ': ' Export{k,end}]);
            elseif length(ParamFN) == 2
                fprintf(fid,'\r\n%s',[Export{k,1} '_' Export{k,2} '_' Export{k,3} '_' ...
                    Export{k,4} ': ' Export{k,end}]);
            elseif length(ParamFN) == 3
                fprintf(fid,'\r\n%s',[Export{k,1} '_' Export{k,2} '_' Export{k,3} '_' ...
                    Export{k,4} '_' Export{k,5} ': ' Export{k,end}]);
            end
        end
    end
end

% Closing the file
fclose(fid);