function [modelGutUniversal] = create_model_from_file(filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define host-gut model with contributions of host and bacteria
% to drug metabolism


%% Import data from text file.
delimiter = ';';
% check the number of columns fo further import
fid = fopen(filename,'rt');
tLines = fgets(fid);
numCols = numel(strfind(tLines,delimiter)) + 1;
fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = [repmat('%s', 1,numCols) '%[^\n\r]'];
%% Open the text file.
fileID = fopen(filename,'r');
%% Read columns of data according to format string.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
%% Close the text file.
fclose(fileID);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
search_for_columns = {'Model entity', 'Name', 'Description', 'Units',...
                      'Initial value', 'Fit data'};
search_for_columns_idx = zeros(size(search_for_columns));
header_idx = zeros(size(search_for_columns));
for i=1:numCols
    find_header = find(cellfun(@(x) ismember(x,search_for_columns), dataArray{i}));
    if length(find_header) > 1 %more than one word in column corresponds to column names
        disp(sprintf(['Error: more than one column keyword\n[' ...
                      repmat('%s | ', 1,length(search_for_columns)-1) ...
                      '%s]\nin column %d\n'],...
                      search_for_columns{:}, i));
    end
    if length(find_header)==1
        header_idx(i) = find_header;
        search_for_columns_idx(i) = find(ismember(search_for_columns, dataArray{i}(find_header)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check that all the headers are on one line
if length(unique(header_idx))>1
    disp(sprintf(['Error: headers are not on the same line\n['...
                   repmat('%s | ', 1,length(search_for_columns)-1) ...
                   '%s]\n[%d | %d | %d | %d | %d]\n'],...
            search_for_columns{:}, header_idx(search_for_columns_idx)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check that all the headers are present
if nnz(search_for_columns_idx) ~= length(search_for_columns_idx)
    missing_columns = repmat('%s\n', 1,nnz(search_for_columns_idx==0));
    disp(sprintf(['Error: missing column(s):\n' missing_columns '\n'],...
            search_for_columns{search_for_columns_idx==0}));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove lines before headers
for i=1:length(search_for_columns_idx)
    dataArray{search_for_columns_idx(i)} = ...
        dataArray{search_for_columns_idx(i)}(header_idx(search_for_columns_idx(i)):end);
end
% exchange empty with ''
dataArray{search_for_columns_idx(6)}(cellfun(@(x) isempty(x),...
            dataArray{search_for_columns_idx(6)})) = {''};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse model species, reactions, coefficients and constants
modelEntities = dataArray{search_for_columns_idx(1)};
modelReactions = find(ismember(modelEntities, 'Equation'));
modelSpecies = find(ismember(modelEntities, 'Metabolite'));
modelParameters = find(ismember(modelEntities, 'Parameter'));
modelConstants = find(ismember(modelEntities, 'Constant'));
modelUnknown = setdiff(modelEntities, {'Equation', 'Metabolite', 'Parameter', 'Constant'});
modelUnknown = setdiff(modelUnknown, search_for_columns);
if ~isempty(modelUnknown)
    disp(sprintf(['Warning: unknown model entity:\n', repmat('%s\n', 1,length(modelUnknown))],...
         modelUnknown{:}));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the model
% define model name as file name (remove extension and folder)
modelName = strsplit(filename, {filesep, '.'});
modelName = modelName{end-1};
modelGutUniversal = sbiomodel(modelName);

%% Add Species and Set initial concentrations of species
for i=1:length(modelSpecies)
    addspecies(modelGutUniversal, dataArray{search_for_columns_idx(2)}{modelSpecies(i)},...
               'InitialAmount',str2double(dataArray{search_for_columns_idx(5)}{modelSpecies(i)}),...
               'InitialAmountUnit', dataArray{search_for_columns_idx(4)}{modelSpecies(i)});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add reations     
for i=1:length(modelReactions)
    addreaction(modelGutUniversal, dataArray{search_for_columns_idx(2)}{modelReactions(i)},...
        'ReactionRate', dataArray{search_for_columns_idx(3)}{modelReactions(i)});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add parameters     
for i=1:length(modelParameters)
   % define the parameter value
    addparameter(modelGutUniversal, dataArray{search_for_columns_idx(2)}{modelParameters(i)},...
                 str2double(dataArray{search_for_columns_idx(5)}{modelParameters(i)}),...
                 'ValueUnits', dataArray{search_for_columns_idx(4)}{modelParameters(i)},...
                 'Notes', dataArray{search_for_columns_idx(6)}{modelParameters(i)});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add constants as parameters     
for i=1:length(modelConstants)
   % define the parameter value
    addparameter(modelGutUniversal, dataArray{search_for_columns_idx(2)}{modelConstants(i)},...
                 str2double(dataArray{search_for_columns_idx(5)}{modelConstants(i)}),...
                 'ValueUnits', dataArray{search_for_columns_idx(4)}{modelConstants(i)},...
                 'Notes', dataArray{search_for_columns_idx(6)}{modelConstants(i)},...
                 'Tag', 'constant');
end