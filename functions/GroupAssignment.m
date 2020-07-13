function [Output,GroupAssign]=GroupAssignment(Data,BetweenLevels,WithinLevels,varargin)
% GroupAssignment() - Automatically assign data to their respective groups/
%                 conditions based on their names.
% Usage:
%          >> [Output,GroupAssign] = GroupAssignment( data, 'key','val'... );
%
% Inputs:
%   Data          = XX.
%   BetweenLevels = e.g. {'Group1, Group2'} XX {default: {''} 
%   WithinLevels  = e.g. {'Condition1, Condition2'} XX {default: {''} 

% Optional inputs:

% Outputs:
%   Output        = XX
%   GroupAssign   = XX will be empty if no between-subject factors.

%% Set defaults
nargs = nargin;
if nargs > 3
    for i = 1:2:length(varargin)
        Param = lower(varargin{i});
        Value = varargin{i+1};
        switch Param
            case 'numcontour'
                CONTOURNUM = Value;
            case 'electrodes'
                ELECTRODES = lower(Value);
        end
    end
end


ListSubj = fieldnames(Data);
Idx=ones(1,length(BetweenLevels)*length(WithinLevels));
GroupAssign=zeros(length(ListSubj),1);
Order = [repmat(1:length(BetweenLevels),[1,length(WithinLevels)])' ...
    reshape(repmat(1:length(WithinLevels),[length(BetweenLevels),1]),...
    [length(BetweenLevels)*length(WithinLevels),1])];

%% ASSIGNMENTS 

% ONLY BETWEEN OR WITHIN-BETWEEN DESIGNS
% Between-subject (BS)
if  exist('BetweenLevels','var') 
    
    for s=1:length(ListSubj)
        
        % For each BS levels
        for k=1:length(BetweenLevels)
            
            if contains(ListSubj{s},BetweenLevels{k}) 
                
                % Position for index
                Pos = k;
                
                % Within-subject (WS)
                if exist('WithinLevels','var') 
                    
                    % For each BS levels
                    for m=1:length(WithinLevels)
                        
                        if contains(ListSubj{s},WithinLevels{m}) 
                            
                            % Retrieving position
                            Pos = find(Order(:,1)==k & Order(:,2)==m);
                            
                            % Fill the database
                            Output.(sprintf('Group%d',k)).(sprintf('Condition%d',m))(:,:,Idx(Pos))=Data.(ListSubj{s});

                            % Update the index
                            Idx(Pos)=Idx(Pos)+1;
                        end

                    end
                    
                    % Prepare the group assignement array
                    GroupAssign(s,1)=k;
                    
                else
                    % Fill the database
                    Output.(sprintf('Group%d',k))(:,:,Idx(k))=Data.(ListSubj{s});

                    % Update the index
                    Idx(Pos)=Idx(Pos)+1;

                    % Prepare the group assignement array
                    GroupAssign(s,1)=1;
                end
            end
        end
    end
    
% ONLY WITHIN DESIGNS
else
    for s=1:length(ListSubj)
        
        % For each BS levels
        for m=1:length(WithinLevels)
            
            % Position for index
            Pos = m;

            if contains(ListSubj{s},WithinLevels{m}) 

                % Fill the database
                Output.(sprintf('Condition%d',m))(:,:,Idx(Pos))=Data.(ListSubj{s});

                % Update the index
                Idx(Pos)=Idx(Pos)+1;
            end
        end
    end  
end

end