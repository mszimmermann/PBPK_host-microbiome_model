function [t,metNamesMap, gd, useForFitting] = load_data_from_file(filename)
%% Import data from text file.
%% Initialize variables.
delimiter = ';';
startRow = 2;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r','n','UTF-8');
% Skip the BOM (Byte Order Mark).
fseek(fileID, 3, 'bof');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,4,5]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end



%% Create output variable
exampledata = raw;
%% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R;

exampledata(:,1) = cellfun(@(x) x{1}, exampledata(:,1), 'unif', 0);
exampledata(cellfun(@(x) isequal(x, 'NaN'), exampledata(:,4)),4) = {NaN};
exampledata(cellfun(@(x) isequal(x, 'NaN'), exampledata(:,5)),5) = {NaN};

exampledata_species = unique(exampledata(:,1));

% get the time from column 3
exampledata_time = cellfun(@(x) str2double(x), exampledata(:,3));
exampledata_time_unique = unique(exampledata_time);

% combine data into a matrix with time versus species
curDataMatrix = zeros(length(exampledata_time_unique), length(exampledata_species));
useForFitting = zeros(length(exampledata_species),1); 
% convert data from concebtrations into amounts given the volumes in column 2
for i=1:length(exampledata_species)
    curData = cellfun(@(x) str2double(x),exampledata(ismember(exampledata(:,1),exampledata_species{i}),4)) .* ...
              cellfun(@(x) str2double(x),exampledata(ismember(exampledata(:,1),exampledata_species{i}),2));
    curTime = cellfun(@(x) str2double(x),exampledata(ismember(exampledata(:,1),exampledata_species{i}),3));
    [~, idxMatrix, idxCur] = intersect(exampledata_time_unique, curTime, 'stable');
    curDataMatrix(idxMatrix,i) = curData(idxCur);
      
    useForFitting(i) = sum(cellfun(@(x) ~isempty(x), exampledata(ismember(exampledata(:,1),exampledata_species{i}),5)))>0;
end
% create a table for symbiology toolbox
t = array2table([exampledata_time_unique...
                curDataMatrix],...
                'VariableNames', [{'Time'},...
                exampledata_species']);
metNamesMap = [exampledata_species exampledata_species];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert data to SimBiology format
gd = groupedData(t);
useForFitting = find(useForFitting);
