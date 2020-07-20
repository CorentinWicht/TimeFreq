%% FIGURES

% Version 0.4 / 13.07.2020 

%% ------------ XXX ---------- %%

%% Author(s)

% Corentin Wicht (script, protocol) 
% GitHub : https://github.com/CorentinWicht

% If you have questions or want to contribute to this pipeline, feel free 
% to contact corentin.wicht@unifr.ch

% This work is licensed under a Creative Commons Attribution-NonCommercial
% 4.0 International License (CC BY-NC)

clear variables 
%% This script enable to repeat the plotting operation while adjusting parameters

% TO DO:

%- Implement for ANOVAs and other stats designs ! 
% + one-sided tests


addpath(genpath(strcat(pwd,'/Functions/Functions')));
addpath(strcat(pwd,'/Functions/eeglab14_1_2b'));

% Loading the excel output file
[file,path] = uigetfile([pwd '/.xlsx'],'Select the excel file to load');
Table = readtable([path file]);

% Loading the statistics RawTFData (mat files)
path = uigetdir(pwd,'Select the folder containing STATS data (Stats folder)');
AllFiles = dir ([path '\**/*' '.mat']);
StatsFiles = AllFiles(contains({AllFiles.folder},'Stats'));
for k=1:length(StatsFiles)
    TempNames = strsplit(StatsFiles(k).name,'.');
    StatsData.(TempNames{1}) = load([StatsFiles(k).folder '\' StatsFiles(k).name]);
end
StatsDataFN = fieldnames(StatsData);

% Loading the raw data files (mat files)
[file,path] = uigetfile([pwd '*.mat'],'Select the file containing RAW data (RawTFData.mat in Raw folder)');
RawTFData = load([path file]);

% Loading the parameters
ParamFile = dir (['**/*' 'TimeFreq_Param.mat']);
load([ParamFile.folder '\' ParamFile.name])

% Table headers
TableHeaders = Table.Properties.VariableNames;

% Retrieve data from the table
Analyses = table2cell(Table(:,contains(TableHeaders,'Analysis')));
FreqBands = table2cell(Table(:,contains(TableHeaders,'FrequencyBand')));
Level1 = table2cell(Table(:,contains(TableHeaders,'Level1')));
Level2 = table2cell(Table(:,contains(TableHeaders,'Level2')));
Level3 = table2cell(Table(:,contains(TableHeaders,'Level3')));

% For t-tests
if contains(StatsData.(StatsDataFN{1}).Stats.Test,'t-test')
    Idx = contains(TableHeaders,'ClusterThreshold');
    if sum(Idx) == 1
        Threshold = table2array(Table(:,Idx));
        SigClust = table2array(Table(:,contains(TableHeaders,'SigClusters')));
    end
    
% For ANOVAs
elseif contains(StatsData.(StatsDataFN{1}).Stats.Test,'ANOVA')
    Threshold = table2array(Table(:,contains(TableHeaders,'ClusterThreshold')));
    SigClust = table2array(Table(:,contains(TableHeaders,'SigClusters')));
end

% Statistics specifications
Pval = str2double(PromptInputs{5}); 
Normalization = upper(answer{8});

% Frequency bands
BandsList=BandsList(all(~cellfun('isempty', BandsList),2),:);

% Creating the Design structure separating Within and Between factors
Design=[];
for k=1:size(DesignList,2) % (MAY CHANGE !!!)
    if strcmpi(DesignList{1,k},'W')
        Design.Within=DesignList(~cellfun('isempty',DesignList(:,k)),k);
    elseif strcmpi(DesignList{1,k},'B')
        Design.Between=DesignList(~cellfun('isempty',DesignList(:,k)),k);
    end
end

% % Removing empty lines
% DesignFN=fieldnames(Design);
% for k=1:length(DesignFN)
%     FactorsFN = fieldnames(Design.(DesignFN{k}));
%     for j=1:length(FactorsFN)
%         Design.(DesignFN{k}).(FactorsFN{j})=Design.(DesignFN{k}).(FactorsFN{j}) ...
%             (~cellfun('isempty',Design.(DesignFN{k}).(FactorsFN{j})));
%     end
% end

% Starting EEGLAB to retrieve one example file
eeglab 
close gcf

% List of EEG files
EEGFiles = dir([data_folder '\**/*' '.set']);

% Loading the first file as example
EEG = pop_loadset('filename',EEGFiles(1).name,'filepath',EEGFiles(1).folder);

% Average referencing Cz
% if strcmpi(answer{6},'Y')
%     EEG = average_ref(EEG,EEG.chaninfo.nodatchans);
% end

% Retrieving labels
Labels = {EEG.chanlocs.labels};

%% Looping over each analysis

% Until user is satisfied
FirstRun = 1;
DispText = '';
Num = 1;

for k=1:size(Table,1) 
             
    % Retrieving the condition labels
    if isfield(Design,'Between')
        BetweenLevels = Design.Between(3:end);
    end
    if isfield(Design,'Within')
        WithinLevels = Design.Within(3:end);
    end
    CurrentBand = FreqBands{k};

    % Listing all subjects + loading the frequency bins
    TempData = RawTFData.(Analyses{k}).(Level1{k}).(Level2{k}).Data;
    SubjFN = fieldnames(TempData);
    
    AllStats = StatsData.([Analyses{k} '_' Level1{k} '_' Level2{k} '_' FreqBands{k}]).Stats.TFCE;  
    AllFreqs = RawTFData.(Analyses{k}).(Level1{k}).(Level2{k}).Freqs.(SubjFN{1});

    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1) T-TESTS
    
    if contains(StatsData.(StatsDataFN{1}).Stats.Test,'t-test')
    
        % Thresholding the statistics based on alpha (e.g. 0.05)
        AllStats.P_Values(AllStats.P_Values>Pval) = 0;
        AllStats.Obs(AllStats.P_Values==0) = 0;

        % Only if significant clusters survived thresholding
        if SigClust(k)>0 

            % Looking for index of frequency boundaries
            CurrentBounds=str2num(BandsList{contains(BandsList(:,1),CurrentBand),2});
            BoundsIdx=find(ge(AllFreqs,CurrentBounds(1)) & le(AllFreqs,CurrentBounds(2))); 

            % Assigning each participant data to its condition
            p=[1 1];
            for m=1:length(SubjFN)
                if contains(SubjFN{m},BetweenLevels{1})
                    GroupRawData.(BetweenLevels{1})(p(1),:,:) = squeeze(mean(TempData.(SubjFN{m})(:,BoundsIdx,:),2));
                    p(1)=p(1)+1;
                elseif contains(SubjFN{m},BetweenLevels{2})
                    GroupRawData.(BetweenLevels{2})(p(2),:,:)= squeeze(mean(TempData.(SubjFN{m})(:,BoundsIdx,:),2));  
                    p(2)=p(2)+1;
                end
            end

            % Averaging over participants
            GroupRawData = structfun(@(x) squeeze(mean(x,1)),GroupRawData,'UniformOutput',false);
            %% Settings

            % SHOULD LET THE POSSIBILITY TO DEFINE THE RANGE FOR EACH TOPOPLOT.
            % Hence, at the end of the loop, if re-run, clear the variables
            % link the time range and below only create them is doesn't exit
            % (hence, if prompt here is empty for time range).

            % Prompt allowing to adjust on the first figure the colormap
            % options
            ColorMap={};
            % Only for the 1st figure
            while FirstRun
                PromptSetup = {['1) Enter the name of the ColorMap to use for MEANS:'...
                    newline 'e.g. parula, jet, hsv, hot, bluewhitered'...
                    newline 'If left empty, opens a software to design it manually!'],...
                    ['2) Enter the name of the ColorMap to use for STATS and TOPOPLOTS:'...
                    newline 'If left empty, opens a software to design it manually!']...
                    'Enter the number of Color bins:'...
                    'Range of topoplots (is ms, optional)'};
                PromptInputs = inputdlg(PromptSetup,'Figure ColorMap Settings',1,{'jet','hot','50',''});


                % Manual configuration of the color maps
                % COLORMAP 1
                if isempty(PromptInputs{1}) && isempty(ColorMap) || ...
                        isempty(PromptInputs{1}) && strcmpi(ColorMap,'Second ColorMap')
                    if exist('ColParam1','var')
                        [ColMap1,ColParam1] = cubehelix_view(size(ColMap1,1),ColParam1(1),...
                            ColParam1(2),ColParam1(3),ColParam1(4),ColParam1(5:6),ColParam1(7:8));
                    else
                        [ColMap1,ColParam1] = cubehelix_view(str2double(PromptInputs{3}));
                    end
                    % Reversing the direction of the Colormap
                    ColMap1 = flipud(ColMap1);
                end 
                % COLORMAP 2
                if isempty(PromptInputs{2}) && isempty(ColorMap) || ...
                        isempty(PromptInputs{1}) && strcmpi(ColorMap,'First ColorMap')
                    if exist('ColParam2','var')
                        [ColMap2,ColParam2] = cubehelix_view(size(ColMap2,1),ColParam2(1),...
                            ColParam2(2),ColParam2(3),ColParam2(4),ColParam2(5:6),ColParam2(7:8));
                    else
                        [ColMap2,ColParam2] = cubehelix_view(str2double(PromptInputs{3}));
                    end
                    % Reversing the direction of the Colormap
                    ColMap2 = flipud(ColMap2);
                end 
                
                % Determine if dependent/independent samples t-test
                if exist('BetweenLevels','var')
                    CurrentLevels = BetweenLevels;
                elseif exist('WithinLevels','var')
                    CurrentLevels = WithinLevels;
                end

                %% PLOTTING  

                % plot tresholded results
                Fig = figure; clf

                % More space if no bar
                set(gcf, 'MenuBar', 'None', 'ToolBar','None')

                % Grand title over all subplots
                Str1=sprintf('%s-%s-%s: %s [%d-%d Hz]',Analyses{k},Level1{k},Level2{k},...
                    FreqBands{k},CurrentBounds(1),CurrentBounds(end));
                Str2=sprintf('%s_%s_%s_%s',Analyses{k},Level1{k},Level2{k},FreqBands{k});
                annotation('textbox',[.35 .5 .5 .5],'String',...
                Str1,'FitBoxToText','on','FontSize',20,'EdgeColor','w');

                % Data to plot
                DataToPlot1=GroupRawData.(CurrentLevels{1});
                DataToPlot2=GroupRawData.(CurrentLevels{2});
                MinVal=min([min(min(min(DataToPlot1))) min(min(min(DataToPlot2)))]);
                MaxVal=max([max(max(max(DataToPlot1))) max(max(max(DataToPlot2)))]);

                % Axes
                XTStr = round(EEG.xmin*1000):50:round(EEG.xmax*1000,-1);
                Ticks = linspace(1,size(DataToPlot1,2),size(XTStr,2));

                % ColorMap
                if ~exist('ColMap1','var')
                    ColMap1 = [PromptInputs{1} '(' PromptInputs{3} ')'];
                end
                if ~exist('ColMap2','var')
                    ColMap2 = [PromptInputs{2} '(' PromptInputs{3} ')'];
                end

                % 1) Original power of group 1 (mean)
                subplot(4,4,[1 2 5 6])
                imagesc(DataToPlot1)
                xlabel('Time (ms)'), ylabel('Channels')
                hold on; contour(logical(AllStats.Obs),1,'linecolor','k'); hold off;
                title(sprintf('Power for group: %s',strrep(CurrentLevels{1},'_','-'))) 
                set(gca,'ydir','norm')
                set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
                set(gca,'YTick',1:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(1:2:end));

                % COLORBAR
                colormap(gca,ColMap1)
                hcb1=colorbar;
                if strcmpi(Normalization,'Y')
                    set(get(hcb1,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                else
                    set(get(hcb1,'XLabel'),'String','Power (\muV)');                
                end
                caxis([MinVal MaxVal])
                grid on

                 % 2) Original power of group 2 (mean)
                subplot(4,4,[3 4 7 8])
                imagesc(DataToPlot2)
                xlabel('Time (ms)'), ylabel('Channels')
                hold on; contour(logical(AllStats.Obs),1,'linecolor','k'); hold off;
                title(sprintf('Power for group: %s',strrep(CurrentLevels{2},'_','-')))
                set(gca,'ydir','norm')
                set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
                set(gca,'YTick',2:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(2:2:end));

                % COLORBAR
                colormap(gca,ColMap1)
                hcb2=colorbar;
                caxis([MinVal MaxVal])
                if strcmpi(Normalization,'Y')
                    set(get(hcb2,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                else
                    set(get(hcb2,'XLabel'),'String','Power (\muV)');                
                end
                grid on

                % 3) Permutation-corrected t-values
                subplot(4,4,[9 10 13 14])
                h=imagesc(AllStats.Obs);
                xlabel('Time (ms)'), ylabel('Channels')
                title(sprintf('TFCE permutation-corrected test [threshold for sig : %.2f]',...
                    Threshold(k)))
                set(gca,'ydir','normal') 
                set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
                set(gca,'YTick',1:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(1:2:end));

                % COLORBAR
                % Set 0 (i.e. non-sig) values to transparent : 
                set(h,'alphadata',AllStats.Obs~=0)
                cmap = colormap(gca,ColMap2);
                hcb3=colorbar;

                % Everything below/above the t-critical should be white
                % BETTER METHOD BUT UNSURE WILL WORK WITH 2-SIDED EFFECTS !!!!!
                StepSize = (hcb3.Limits(2)-hcb3.Limits(1))/size(cmap,1);
                Steps = hcb3.Limits(2):-StepSize:hcb3.Limits(1);
                Idx=Steps<abs(Threshold(k)) & Steps>-abs(Threshold(k));
                % cmap(flip(Idx(2:end)),:)=repmat([1 1 1],[sum(Idx)-1,1]); % WHY -1 ?? Works if no -1 (for negative and positive values)
                caxis([min(Steps(~Idx)) max(Steps(~Idx))+500])
                % I ADDED THE + 500 otherwise couldn't see the results on top!!!!
                colormap(gca,cmap);
                title(hcb3, 'TFCE-values');
                grid on

                % Timing of significant electrodes
                PosTopoplot = [11 12 15 16];
                [SigElect,SigPoints]=find(AllStats.P_Values~=0);
                SigPoints = unique(SigPoints);
                SigElect = unique(SigElect);
                SigLength = (SigPoints(end) - SigPoints(1))/2;
                ExactPoints = EEG.xmin*EEG.srate:EEG.xmax*EEG.srate;
                ExactTimes = round(SigPoints(1):SigLength:SigPoints(end));

                % Retrieve the most significant electrode/time point
                [~, MostSigPoints] = find(AllStats.P_Values==min(min(AllStats.P_Values(AllStats.P_Values~=0))));    

                % Min/max values for the color bar
                MinVal(1) = min(min(DataToPlot1(:,MostSigPoints(1):MostSigPoints(end))));
                MinVal(2) = min(min(DataToPlot2(:,MostSigPoints(1):MostSigPoints(end))));
                MaxVal(1) = max(max(DataToPlot1(:,MostSigPoints(1):MostSigPoints(end))));
                MaxVal(2) = max(max(DataToPlot2(:,MostSigPoints(1):MostSigPoints(end))));

                % for each topoplot
                for f=1:4
                    % 4) Prepare the subplots
                    subplot(4,4,PosTopoplot(f))

                    % Time in titles
                    if isempty(PromptInputs{4})
                        MinTime = round(ExactPoints(MostSigPoints(1))/EEG.srate*1000);
                        MaxTime = round(ExactPoints(MostSigPoints(end))/EEG.srate*1000,-1);
                    else
                        Temp = str2num(PromptInputs{4});
                        MinTime = Temp(1);
                        MaxTime = Temp(end);
                        MostSigPoints = ExactPoints(round(MinTime/EEG.srate*1000)): ...
                            ExactPoints(round(MaxTime/EEG.srate*1000)); % NOT SURE ABOUT THIS !!! 
                    end

                    if f==1
                        % Topoplots
                        topoplotIndie(DataToPlot1(:,MostSigPoints(1):MostSigPoints(end)),...
                            EEG.chanlocs,'sigelect',SigElect,'colsigelect','k');
                        title(sprintf('%s [%d-%d ms]',strrep(CurrentLevels{1},'_','-'),...
                            MinTime,MaxTime))

                        % Colorbar
                        cmap = colormap(gca,ColMap1);
                        hcb4=colorbar;
                        caxis([min(MinVal) max(MaxVal)])
                        if strcmpi(Normalization,'Y')
                            set(get(hcb4,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                        else
                            set(get(hcb4,'XLabel'),'String','Power (\muV)');                
                        end
                    elseif f==2
                        % Topoplots
                        topoplotIndie(DataToPlot2(:,MostSigPoints(1):MostSigPoints(end)),...
                            EEG.chanlocs,'sigelect',SigElect,'colsigelect','k');
                        title(sprintf('%s [%d-%d ms]',strrep(CurrentLevels{2},'_','-'),...
                            MinTime,MaxTime))

                        % Colorbar
                        cmap = colormap(gca,ColMap1);
                        hcb4=colorbar;
                        caxis([min(MinVal) max(MaxVal)])
                        if strcmpi(Normalization,'Y')
                            set(get(hcb4,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                        else
                            set(get(hcb4,'XLabel'),'String','Power (\muV)');
                        end
                    else
                        % Topoplots
                        topoplotIndie(DataToPlot1(:,ExactTimes(f-2):ExactTimes(f-1))- ...
                            DataToPlot2(:,ExactTimes(f-2):ExactTimes(f-1)),EEG.chanlocs,'sigelect',SigElect,'colsigelect','k');
                        title(sprintf('\\Delta %s - %s  [%d-%d ms]',...
                            strrep(CurrentLevels{1},'_','-'),strrep(CurrentLevels{2},'_','-'),...
                            ExactPoints(round(ExactTimes(f-2)/EEG.srate*1000)), ...
                            ExactPoints(round(ExactTimes(f-1)/EEG.srate*1000))))

                        % Colorbar
                        cmap = colormap(gca,ColMap2);
                        hcb4=colorbar;
                        if strcmpi(Normalization,'Y')
                            set(get(hcb4,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                        else
                            set(get(hcb4,'XLabel'),'String','Power (\muV)');
                        end
                    end
                end

                % Set figure as fullscreen
                set(Fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

                %% Export 

                % Exporting figures

                % Export input
                Export = questdlg('What would you like to do?', 'Figure export',...
                    'Export in .fig', 'Export in .bmp','Do not export','Export in .bmp');

                if ~strcmpi(Export,'Do not export')
                    % Create export folders
                    mkdir ([save_folder '/Figures'])
                end

                % Are you satisfied?
                Satisfy = questdlg('Are you satisfied with your colormap choice & figure export format specifications ?', 'Processing all figures',...
                    'yes', 'no','no');
                if strcmpi(Satisfy,'no') && exist('ColParam1','var') || ...
                        strcmpi(Satisfy,'no') && exist('ColParam2','var')
                    ColorMap = questdlg(['Would you like to keep the ColorMap that you defined manually ?'...
                        newline 'IF NOT, CLOSE THE WINDOW'], 'Manual ColorMaps',...
                        'First ColorMap', 'Second ColorMap','Both ColorMaps','Both ColorMaps');
                end
                if strcmpi(Satisfy,'yes')
                    FirstRun = 0;

                    % Exporting 
                    if strcmpi(Export,'Export in .bmp')
                        SaveFigures(gcf,[save_folder '/Figures/' Str2],'w','bmp');
                    elseif strcmpi(Export,'Export in .fig')
                        savefig([save_folder '/Figures/' Str2])
                        close gcf
                    elseif strcmpi(Export,'Do not export')
                        close gcf
                    end

                    % Free space
                    clear GroupRawData
                elseif strcmpi(Satisfy,'no')
                    close gcf
                elseif isempty(Satisfy)
                    disp('The program has been stopped'); 
                    break
                end
            end
        end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2) ANOVAs ( CURRENTLY ONLY FOR MIXED 2x2 ! )
    elseif contains(StatsData.(StatsDataFN{1}).Stats.Test,'Mixed ANOVA')
        
        % For both main effects + interaction
        % A = Within-subject main
        % B = Between-subject main
        % AB = Interaction
        EffectsStr = fieldnames(AllStats);
        for t=1:3 
            
            % Current statistic
            CurrentStats = AllStats.(EffectsStr{t});
        
            % Thresholding the statistics based on alpha (e.g. 0.05)
            CurrentStats.P_Values(CurrentStats.P_Values>Pval) = 0;
            CurrentStats.Obs(CurrentStats.P_Values==0) = 0;

            % Only if significant clusters survived thresholding
            if SigClust(k,t)>0 

                % Looking for index of frequency boundaries
                CurrentBounds=str2num(BandsList{contains(BandsList(:,1),CurrentBand),2});
                BoundsIdx=find(ge(AllFreqs,CurrentBounds(1)) & le(AllFreqs,CurrentBounds(2))); 

                % Assigning each participant data to its group/condition
                % Columns are Between-subject factors (B)
                % Lines are Within-subject (A)
                TempDataCell = struct2cell(TempData);
                for m=1:length(BetweenLevels)
                   for l=1:length(WithinLevels)
                        Idx = find(contains(SubjFN,BetweenLevels{m})+contains(SubjFN,WithinLevels{l}) == 2);
                        DataTemp = cellfun(@(x) squeeze(mean(x(:,BoundsIdx,:),2)),TempDataCell(Idx),'UniformOutput',0);
                        for u=1:length(DataTemp)
                            GroupRawData{m,l}(u,:,:) = DataTemp{u};
                        end
                   end
                end
                
                % Averaging over participants
                GroupRawData = cellfun(@(x) squeeze(mean(x,1)),GroupRawData,'UniformOutput',false);
                %% Settings

                % SHOULD LET THE POSSIBILITY TO DEFINE THE RANGE FOR EACH TOPOPLOT.
                % Hence, at the end of the loop, if re-run, clear the variables
                % link the time range and below only create them is doesn't exit
                % (hence, if prompt here is empty for time range).

                % Prompt allowing to adjust on the first figure the colormap
                % options
                ColorMap={};
                % Only for the 1st figure
                while FirstRun
                    PromptSetup = {['1) Enter the name of the ColorMap to use for MEANS:'...
                        newline 'e.g. parula, jet, hsv, hot, bluewhitered'...
                        newline 'If left empty, opens a software to design it manually!'],...
                        ['2) Enter the name of the ColorMap to use for STATS and TOPOPLOTS:'...
                        newline 'If left empty, opens a software to design it manually!']...
                        'Enter the number of Color bins:'...
                        'Range of topoplots (is ms, optional)'};
                    PromptInputs = inputdlg(PromptSetup,'Figure ColorMap Settings',1,{'jet','hot','50',''});


                    % Manual configuration of the color maps
                    % COLORMAP 1
                    if isempty(PromptInputs{1}) && isempty(ColorMap) || ...
                            isempty(PromptInputs{1}) && strcmpi(ColorMap,'Second ColorMap')
                        if exist('ColParam1','var')
                            [ColMap1,ColParam1] = cubehelix_view(size(ColMap1,1),ColParam1(1),...
                                ColParam1(2),ColParam1(3),ColParam1(4),ColParam1(5:6),ColParam1(7:8));
                        else
                            [ColMap1,ColParam1] = cubehelix_view(str2double(PromptInputs{3}));
                        end
                        % Reversing the direction of the Colormap
                        ColMap1 = flipud(ColMap1);
                    end 
                    % COLORMAP 2
                    if isempty(PromptInputs{2}) && isempty(ColorMap) || ...
                            isempty(PromptInputs{1}) && strcmpi(ColorMap,'First ColorMap')
                        if exist('ColParam2','var')
                            [ColMap2,ColParam2] = cubehelix_view(size(ColMap2,1),ColParam2(1),...
                                ColParam2(2),ColParam2(3),ColParam2(4),ColParam2(5:6),ColParam2(7:8));
                        else
                            [ColMap2,ColParam2] = cubehelix_view(str2double(PromptInputs{3}));
                        end
                        % Reversing the direction of the Colormap
                        ColMap2 = flipud(ColMap2);
                    end 

                    %% PLOTTING  

                    % plot tresholded results
                    Fig = figure; clf

                    % More space if no bar
                    set(gcf, 'MenuBar', 'None', 'ToolBar','None')

                    % Grand title over all subplots
                    Str1=sprintf('%s-%s-%s-%s: %s [%d-%d Hz]',Analyses{k},Level1{k},Level2{k},...
                        ['Contrast' EffectsStr{t}],FreqBands{k},CurrentBounds(1),CurrentBounds(end));
                    Str2=sprintf('%s_%s_%s_%s_%s',Analyses{k},Level1{k},Level2{k},['Contrast' EffectsStr{t}],FreqBands{k});
                    annotation('textbox',[.35 .5 .5 .5],'String',...
                    Str1,'FitBoxToText','on','FontSize',20,'EdgeColor','w');

                    % DATA TO PLOT + SETTINGS
                    % Main Effects
                    if t == 1 % Within
                        DataToPlot1=mean(reshape([GroupRawData{1,1} GroupRawData{1,2}],...
                            [size(GroupRawData{1,1},1) size(GroupRawData{1,1},2) 2]),3);
                        DataToPlot2=mean(reshape([GroupRawData{2,1} GroupRawData{2,2}],...
                            [size(GroupRawData{2,1},1) size(GroupRawData{2,1},2) 2]),3);
                        
                        % Min/max (for scale)
                        MinVal=min([min(min(min(DataToPlot1))) min(min(min(DataToPlot2)))]);
                        MaxVal=max([max(max(max(DataToPlot1))) max(max(max(DataToPlot2)))]);
                        
                        % Define current levels
                        CurrentLevels = WithinLevels;
                        
                        % Position of subplots
                        Plot1 = [1 2 5 6]; Plot2 = [3 4 7 8];
                    elseif t == 2 % Between
                        DataToPlot1=mean(reshape([GroupRawData{1,1} GroupRawData{2,1}],...
                            [size(GroupRawData{1,1},1) size(GroupRawData{1,1},2) 2]),3);
                        DataToPlot2=mean(reshape([GroupRawData{1,2} GroupRawData{2,2}],...
                            [size(GroupRawData{1,2},1) size(GroupRawData{1,2},2) 2]),3);
                        
                        % Min/max (for scale)
                        MinVal=min([min(min(min(DataToPlot1))) min(min(min(DataToPlot2)))]);
                        MaxVal=max([max(max(max(DataToPlot1))) max(max(max(DataToPlot2)))]);
                        
                        % Define current levels
                        CurrentLevels = BetweenLevels;
                        
                        % Number of subplots
                        Plot1 = [1 2 5 6]; Plot2 = [3 4 7 8];
                        
                    % Interaction
                    else 
                        % Within
                        DataToPlot1=mean(reshape([GroupRawData{1,1} GroupRawData{1,2}],...
                            [size(GroupRawData{1,1},1) size(GroupRawData{1,1},2) 2]),3);
                        DataToPlot2=mean(reshape([GroupRawData{2,1} GroupRawData{2,2}],...
                            [size(GroupRawData{2,1},1) size(GroupRawData{2,1},2) 2]),3);
                        
                        % Between
                        DataToPlot3=mean(reshape([GroupRawData{1,1} GroupRawData{2,1}],...
                            [size(GroupRawData{1,1},1) size(GroupRawData{1,1},2) 2]),3);
                        DataToPlot4=mean(reshape([GroupRawData{1,2} GroupRawData{2,2}],...
                            [size(GroupRawData{1,2},1) size(GroupRawData{1,2},2) 2]),3);
                        
                        % Min/max (for scale)
                        MinVal=min([min(min(min(DataToPlot1))) min(min(min(DataToPlot2))) ...
                            min(min(min(DataToPlot3))) min(min(min(DataToPlot4)))]);
                        MaxVal=max([max(max(max(DataToPlot1))) max(max(max(DataToPlot2)))...
                            max(max(max(DataToPlot3))) max(max(max(DataToPlot4)))]);
                        
                        % Number of subplots
                        Plot1 = [1 2]; Plot2 = [3 4];
                        Plot3 = [5 6]; Plot4 = [7 8];
                        
                        % Define current levels (concatenate all)
                        CurrentLevels = [WithinLevels;BetweenLevels];
                    end

                    % Axes
                    XTStr = round(EEG.xmin*1000):50:round(EEG.xmax*1000,-1);
                    Ticks = linspace(1,size(DataToPlot1,2),size(XTStr,2));

                    % ColorMap
                    if ~exist('ColMap1','var')
                        ColMap1 = [PromptInputs{1} '(' PromptInputs{3} ')'];
                    end
                    if ~exist('ColMap2','var')
                        ColMap2 = [PromptInputs{2} '(' PromptInputs{3} ')'];
                    end

                    
                    %% SUBPLOTS
                    StringCompare = {'condition','group'};
                    
                    % for main effects
                    if t<3 
                        %% 1) Observed power of group/condition 1
                        subplot(4,4,Plot1)
                        imagesc(DataToPlot1)
                        xlabel('Time (ms)'), ylabel('Channels')
                        hold on; contour(logical(CurrentStats.Obs),1,'linecolor','k'); hold off;
                        title(sprintf('Power for %s: %s',StringCompare{t},strrep(CurrentLevels{1},'_','-'))) 
                        set(gca,'ydir','norm')
                        set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
                        set(gca,'YTick',1:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(1:2:end));

                        % COLORBAR
                        colormap(gca,ColMap1)
                        hcb1=colorbar;
                        if strcmpi(Normalization,'Y')
                            set(get(hcb1,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                        else
                            set(get(hcb1,'XLabel'),'String','Power (\muV)');                
                        end
                        caxis([MinVal MaxVal])
                        grid on

                        %% 2) Observed power of group/condition 2

                        subplot(4,4,Plot2)
                        imagesc(DataToPlot2)
                        xlabel('Time (ms)'), ylabel('Channels')
                        hold on; contour(logical(CurrentStats.Obs),1,'linecolor','k'); hold off;
                        title(sprintf('Power for %s: %s',StringCompare{t},strrep(CurrentLevels{2},'_','-')))
                        set(gca,'ydir','norm')
                        set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
                        set(gca,'YTick',2:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(2:2:end));

                        % COLORBAR
                        colormap(gca,ColMap1)
                        hcb2=colorbar;
                        caxis([MinVal MaxVal])
                        if strcmpi(Normalization,'Y')
                            set(get(hcb2,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                        else
                            set(get(hcb2,'XLabel'),'String','Power (\muV)');                
                        end
                        grid on
                        
                    % for interaction
                    elseif t == 3
                        %% 1) Observed power for Condition 1
                        subplot(4,4,Plot1)
                        imagesc(DataToPlot1)
                        xlabel('Time (ms)'), ylabel('Channels')
                        hold on; contour(logical(CurrentStats.Obs),1,'linecolor','k'); hold off;
                        title(sprintf('Power for %s: %s',StringCompare{1},strrep(CurrentLevels{1},'_','-'))) 
                        set(gca,'ydir','norm')
                        set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
                        set(gca,'YTick',1:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(1:2:end));

                        % COLORBAR
                        colormap(gca,ColMap1)
                        hcb1=colorbar;
                        if strcmpi(Normalization,'Y')
                            set(get(hcb1,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                        else
                            set(get(hcb1,'XLabel'),'String','Power (\muV)');                
                        end
                        caxis([MinVal MaxVal])
                        grid on

                        %% 2) Observed power for Condition 2

                        subplot(4,4,Plot2)
                        imagesc(DataToPlot2)
                        xlabel('Time (ms)'), ylabel('Channels')
                        hold on; contour(logical(CurrentStats.Obs),1,'linecolor','k'); hold off;
                        title(sprintf('Power for %s: %s',StringCompare{1},strrep(CurrentLevels{2},'_','-')))
                        set(gca,'ydir','norm')
                        set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
                        set(gca,'YTick',2:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(2:2:end));

                        % COLORBAR
                        colormap(gca,ColMap1)
                        hcb2=colorbar;
                        caxis([MinVal MaxVal])
                        if strcmpi(Normalization,'Y')
                            set(get(hcb2,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                        else
                            set(get(hcb2,'XLabel'),'String','Power (\muV)');                
                        end
                        grid on
                        
                        %% 3) Observed power for Group 1
                        subplot(4,4,Plot3)
                        imagesc(DataToPlot3)
                        xlabel('Time (ms)'), ylabel('Channels')
                        hold on; contour(logical(CurrentStats.Obs),1,'linecolor','k'); hold off;
                        title(sprintf('Power for %s: %s',StringCompare{2},strrep(CurrentLevels{3},'_','-'))) 
                        set(gca,'ydir','norm')
                        set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
                        set(gca,'YTick',1:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(1:2:end));

                        % COLORBAR
                        colormap(gca,ColMap1)
                        hcb1=colorbar;
                        if strcmpi(Normalization,'Y')
                            set(get(hcb1,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                        else
                            set(get(hcb1,'XLabel'),'String','Power (\muV)');                
                        end
                        caxis([MinVal MaxVal])
                        grid on

                        %% 4) Observed power for Group 2

                        subplot(4,4,Plot4)
                        imagesc(DataToPlot4)
                        xlabel('Time (ms)'), ylabel('Channels')
                        hold on; contour(logical(CurrentStats.Obs),1,'linecolor','k'); hold off;
                        title(sprintf('Power for %s: %s',StringCompare{2},strrep(CurrentLevels{4},'_','-')))
                        set(gca,'ydir','norm')
                        set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
                        set(gca,'YTick',2:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(2:2:end));

                        % COLORBAR
                        colormap(gca,ColMap1)
                        hcb2=colorbar;
                        caxis([MinVal MaxVal])
                        if strcmpi(Normalization,'Y')
                            set(get(hcb2,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                        else
                            set(get(hcb2,'XLabel'),'String','Power (\muV)');                
                        end
                        grid on
                    end

                    %% 3) Permutation-corrected t-values
                    subplot(4,4,[9 10 13 14])
                    h=imagesc(CurrentStats.Obs);
                    xlabel('Time (ms)'), ylabel('Channels')
                    title(sprintf('TFCE permutation-corrected test [threshold for sig : %.2f]',...
                        Pval))
                    set(gca,'ydir','normal') 
                    set(gca,'XTick',Ticks,'XTickLabel',num2cell(XTStr));
                    set(gca,'YTick',1:2:length({EEG.chanlocs.labels}),'YTickLabel',Labels(1:2:end));

                    % COLORBAR
                    % Set 0 (i.e. non-sig) values to transparent : 
                    set(h,'alphadata',CurrentStats.Obs~=0)
                    cmap = colormap(gca,ColMap2);
                    hcb3=colorbar;

                    % Everything below/above the t-critical should be white
                    % BETTER METHOD BUT UNSURE WILL WORK WITH 2-SIDED EFFECTS !!!!!
                    % THIS WILL ONLY WORK IF I REPORT THE TFCE RESULTS...
                    
%                     StepSize = (hcb3.Limits(2)-hcb3.Limits(1))/size(cmap,1);
%                     Steps = hcb3.Limits(2):-StepSize:hcb3.Limits(1);
%                     Idx=Steps<abs(Threshold(k)) & Steps>-abs(Threshold(k));
%                     caxis([min(Steps(~Idx)) max(Steps(~Idx))+500])
                    % I ADDED THE + 500 otherwise couldn't see the results on top!!!!
                    colormap(gca,cmap);
                    title(hcb3, 't-values');
                    grid on

                    % Position to plot the topoplots
                    PosTopoplot = [11 12 15 16];
                    % Timing of significant electrodes
                    [SigElect,SigPoints]=find(CurrentStats.P_Values~=0);
                    SigPoints = unique(SigPoints);
                    SigElect = unique(SigElect);
                    SigLength = (SigPoints(end) - SigPoints(1))/2;
                    ExactPoints = EEG.xmin*EEG.srate:EEG.xmax*EEG.srate;
                    ExactTimes = round(SigPoints(1):SigLength:SigPoints(end));

                    % Retrieve the most significant electrode/time point
                    [~, MostSigPoints] = find(CurrentStats.P_Values==min(min(CurrentStats.P_Values(CurrentStats.P_Values~=0))));    

                    % For main effects
                    if t<3
                        
                        % Min/max values for the color bar
                        MinVal(1) = min(min(DataToPlot1(:,MostSigPoints(1):MostSigPoints(end))));
                        MinVal(2) = min(min(DataToPlot2(:,MostSigPoints(1):MostSigPoints(end))));
                        MaxVal(1) = max(max(DataToPlot1(:,MostSigPoints(1):MostSigPoints(end))));
                        MaxVal(2) = max(max(DataToPlot2(:,MostSigPoints(1):MostSigPoints(end))));

                        % for each topoplot
                        for f=1:4
                            % 4) Prepare the subplots
                            subplot(4,4,PosTopoplot(f))

                            % Time in titles
                            if isempty(PromptInputs{4})
                                MinTime = round(ExactPoints(MostSigPoints(1))/EEG.srate*1000);
                                MaxTime = round(ExactPoints(MostSigPoints(end))/EEG.srate*1000,-1);
                            else
                                Temp = str2num(PromptInputs{4});
                                MinTime = Temp(1);
                                MaxTime = Temp(end);
                                MostSigPoints = ExactPoints(round(MinTime/EEG.srate*1000)): ...
                                    ExactPoints(round(MaxTime/EEG.srate*1000)); % NOT SURE ABOUT THIS !!! 
                            end

                            if f==1
                                % Topoplots
                                topoplotIndie(DataToPlot1(:,MostSigPoints(1):MostSigPoints(end)),...
                                    EEG.chanlocs,'sigelect',SigElect,'colsigelect','k');
                                title(sprintf('%s [%d %dms]',strrep(CurrentLevels{1},'_','-'),...
                                    MinTime,MaxTime))

                                % Colorbar
                                cmap = colormap(gca,ColMap1);
                                hcb4=colorbar;
                                caxis([min(MinVal) max(MaxVal)])
                                if strcmpi(Normalization,'Y')
                                    set(get(hcb4,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                                else
                                    set(get(hcb4,'XLabel'),'String','Power (\muV)');                
                                end
                            elseif f==2
                                % Topoplots
                                topoplotIndie(DataToPlot2(:,MostSigPoints(1):MostSigPoints(end)),...
                                    EEG.chanlocs,'sigelect',SigElect,'colsigelect','k');
                                title(sprintf('%s [%d %dms]',strrep(CurrentLevels{2},'_','-'),...
                                    MinTime,MaxTime))

                                % Colorbar
                                cmap = colormap(gca,ColMap1);
                                hcb4=colorbar;
                                caxis([min(MinVal) max(MaxVal)])
                                if strcmpi(Normalization,'Y')
                                    set(get(hcb4,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                                else
                                    set(get(hcb4,'XLabel'),'String','Power (\muV)');
                                end
                            else
                                % Topoplots
                                topoplotIndie(DataToPlot1(:,ExactTimes(f-2):ExactTimes(f-1))- ...
                                    DataToPlot2(:,ExactTimes(f-2):ExactTimes(f-1)),EEG.chanlocs,'sigelect',SigElect,'colsigelect','k');
                                title(sprintf('\\Delta %s - %s  [%d %dms]',...
                                    strrep(CurrentLevels{1},'_','-'),strrep(CurrentLevels{2},'_','-'),...
                                    ExactPoints(round(ExactTimes(f-2)/EEG.srate*1000)), ...
                                    ExactPoints(round(ExactTimes(f-1)/EEG.srate*1000))))

                                % Colorbar
                                cmap = colormap(gca,ColMap2);
                                hcb4=colorbar;
                                if strcmpi(Normalization,'Y')
                                    set(get(hcb4,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                                else
                                    set(get(hcb4,'XLabel'),'String','Power (\muV)');
                                end
                            end
                        end
                        
                    % For interaction effect                        
                    elseif t == 3
                        
                        % Min/max values for the color bar
                        MinVal(1) = min(min(DataToPlot1(:,MostSigPoints(1):MostSigPoints(end))));
                        MinVal(2) = min(min(DataToPlot2(:,MostSigPoints(1):MostSigPoints(end))));
                        MinVal(3) = min(min(DataToPlot3(:,MostSigPoints(1):MostSigPoints(end))));
                        MinVal(4) = min(min(DataToPlot4(:,MostSigPoints(1):MostSigPoints(end))));
                        
                        MaxVal(1) = max(max(DataToPlot1(:,MostSigPoints(1):MostSigPoints(end))));
                        MaxVal(2) = max(max(DataToPlot2(:,MostSigPoints(1):MostSigPoints(end))));
                        MaxVal(3) = max(max(DataToPlot3(:,MostSigPoints(1):MostSigPoints(end))));
                        MaxVal(4) = max(max(DataToPlot4(:,MostSigPoints(1):MostSigPoints(end))));

                        % for each topoplot
                        for f=1:4
                            % 4) Prepare the subplots
                            subplot(4,4,PosTopoplot(f))

                            % Time in titles
                            if isempty(PromptInputs{4})
                                MinTime = round(ExactPoints(MostSigPoints(1))/EEG.srate*1000);
                                MaxTime = round(ExactPoints(MostSigPoints(end))/EEG.srate*1000,-1);
                            else
                                Temp = str2num(PromptInputs{4});
                                MinTime = Temp(1);
                                MaxTime = Temp(end);
                                MostSigPoints = ExactPoints(round(MinTime/EEG.srate*1000)): ...
                                    ExactPoints(round(MaxTime/EEG.srate*1000)); % NOT SURE ABOUT THIS !!! 
                            end

                            if f<3
                                % Topoplots
                                Pos = 1;
                                topoplotIndie(DataToPlot1(:,ExactTimes(Pos):ExactTimes(Pos+1))- ...
                                    DataToPlot2(:,ExactTimes(Pos):ExactTimes(Pos+1)),EEG.chanlocs,'sigelect',SigElect,'colsigelect','k');
                                title(sprintf('\\Delta %s - %s  [%d %dms]',...
                                    strrep(CurrentLevels{1},'_','-'),strrep(CurrentLevels{2},'_','-'),...
                                    ExactPoints(round(ExactTimes(Pos)/EEG.srate*1000)), ...
                                    ExactPoints(round(ExactTimes(Pos+1)/EEG.srate*1000))))

                                % Colorbar
                                cmap = colormap(gca,ColMap2);
                                hcb4=colorbar;
                                if strcmpi(Normalization,'Y')
                                    set(get(hcb4,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                                else
                                    set(get(hcb4,'XLabel'),'String','Power (\muV)');
                                end
                            else
                                Pos = 2;
                                % Topoplots
                                topoplotIndie(DataToPlot3(:,ExactTimes(Pos):ExactTimes(Pos+1))- ...
                                    DataToPlot4(:,ExactTimes(Pos):ExactTimes(Pos+1)),EEG.chanlocs,'sigelect',SigElect,'colsigelect','k');
                                title(sprintf('\\Delta %s - %s  [%d %dms]',...
                                    strrep(CurrentLevels{3},'_','-'),strrep(CurrentLevels{4},'_','-'),...
                                    ExactPoints(round(ExactTimes(Pos)/EEG.srate*1000)), ...
                                    ExactPoints(round(ExactTimes(Pos+1)/EEG.srate*1000))))

                                % Colorbar
                                cmap = colormap(gca,ColMap2);
                                hcb4=colorbar;
                                if strcmpi(Normalization,'Y')
                                    set(get(hcb4,'XLabel'),'String','Power (10*log10(\muV^2/Hz))');
                                else
                                    set(get(hcb4,'XLabel'),'String','Power (\muV)');
                                end
                            end
                        end
                    end

                    % Set figure as fullscreen
                    set(Fig, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);

                    %% Export 

                    % Exporting figures

                    % Export input
                    Export = questdlg('What would you like to do?', 'Figure export',...
                        'Export in .fig', 'Export in .bmp','Do not export','Export in .bmp');

                    if ~strcmpi(Export,'Do not export')
                        % Create export folders
                        mkdir ([save_folder '/Figures'])
                    end

                    % Are you satisfied?
                    Satisfy = questdlg('Are you satisfied with your colormap choice & figure export format specifications ?', 'Processing all figures',...
                        'yes', 'no','no');
                    if strcmpi(Satisfy,'no') && exist('ColParam1','var') || ...
                            strcmpi(Satisfy,'no') && exist('ColParam2','var')
                        ColorMap = questdlg(['Would you like to keep the ColorMap that you defined manually ?'...
                            newline 'IF NOT, CLOSE THE WINDOW'], 'Manual ColorMaps',...
                            'First ColorMap', 'Second ColorMap','Both ColorMaps','Both ColorMaps');
                    end
                    if strcmpi(Satisfy,'yes')
                        FirstRun = 0;

                        % Exporting 
                        if strcmpi(Export,'Export in .bmp')
                            SaveFigures(gcf,[save_folder '/Figures/' Str2],'w','bmp');
                        elseif strcmpi(Export,'Export in .fig')
                            savefig([save_folder '/Figures/' Str2])
                            close gcf
                        elseif strcmpi(Export,'Do not export')
                            close gcf
                        end

                        % Free space
                        clear GroupRawData
                    elseif strcmpi(Satisfy,'no')
                        close gcf
                    elseif isempty(Satisfy)
                        disp('The program has been stopped'); 
                        break
                    end
                end
            end
        end
    end
end

disp('The code is done !')