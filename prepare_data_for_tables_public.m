%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function prepare_data_for_tables_public(t, curVolumes, modelGutUniversal, resultsHost, ...
    metNamesMap, useForFitting, outputFileName, parInit)
% prepare data for tables
expTime = t.Time;

parameterNames = cell(length(modelGutUniversal.Parameters),1);
for i=1:length(modelGutUniversal.Parameters)
    parameterNames{i} = modelGutUniversal.Parameters(i).Name;
end

configsetObj = getconfigset(modelGutUniversal);
configsetObj.SolverOptions.SensitivityAnalysis = false;

predictedResults = zeros(length(expTime), length(modelGutUniversal.Species));
for i=2:length(expTime)
    set(configsetObj, 'StopTime', expTime(i));
    simData = sbiosimulate(modelGutUniversal);
    predictedResults(i,:) = simData.Data(end,:)./curVolumes;
end
experimentalResults = zeros(size(predictedResults));
for i=1:length(simData.dataNames)
    if nnz(ismember(metNamesMap(:,1), simData.dataNames{i}))
        dataMet = metNamesMap(ismember(metNamesMap(:,1), simData.dataNames{i}),2);
        experimentalResults(:,i) = t{:, dataMet}/curVolumes(i);
    end
end

mseDRUG = zeros(1,size(predictedResults,2));
maeDRUG = zeros(1,size(predictedResults,2));
corrDRUG = zeros(1,size(predictedResults,2));
%figure;
for i=1:size(predictedResults,2)
    mseDRUG(i) = mean(( predictedResults(~isnan(experimentalResults(:,i)),i) - ...
                   experimentalResults(~isnan(experimentalResults(:,i)),i)).^2);
    maeDRUG(i) = mean(abs( predictedResults(~isnan(experimentalResults(:,i)),i) - ...
                          experimentalResults(~isnan(experimentalResults(:,i)),i)));
    corrDRUG(i) = corr(predictedResults(~isnan(experimentalResults(:,i)),i),...
                      experimentalResults(~isnan(experimentalResults(:,i)),i));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % diagnostic plotting of experimental and simulated data
%    subplot(3,5,i)
%     plot(expTime, experimentalResults(:,i))
%     hold on
%     plot(expTime, predictedResults(:,i))
%     title(simData.DataNames{i})
end
% set to nan values where corr is nan
mseDRUG(isnan(corrDRUG)) = nan;
maeDRUG(isnan(corrDRUG)) = nan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save data to file
fid = fopen(outputFileName, 'w');
%fprintf(fid, 'Parameter,Initial value,Model estimate,Standard Error,CI 2.5%%,CI97.5%%\n');
fprintf(fid, 'Parameter,Initial value,Model estimate,Standard Error\n');

for i=1:length(modelGutUniversal.Parameters)
    idx = ismember(resultsHost.ParameterEstimates.Name, parameterNames{i});
    if nnz(idx)
        fprintf(fid, '%s,%.3f,%.3f,%.3f\n',...%.3f,%.3f\n', ...
            modelGutUniversal.Parameters(i).Name,...
            parInit(idx),...
            resultsHost.ParameterEstimates.Estimate(idx),...
            resultsHost.ParameterEstimates.StandardError(idx));%resultsHost.parameterci.Results.ConfidenceInterval(idx,:));
    else
        fprintf(fid, '%s,%.3f,%.3f,%.3f\n',...%.3f,%.3f\n', ...
            modelGutUniversal.Parameters(i).Name,...
            modelGutUniversal.Parameters(i).Value,...
            modelGutUniversal.Parameters(i).Value,...
            0);%,0,0);
    end
end
fprintf(fid, '\n');
fprintf(fid, 'Metabolite,Time,Experimental data,Used for fitting,Model prediction,MAE,MSE,PCC\n');
for i=1:size(predictedResults,2)
    if ismember(simData.dataNames{i}, metNamesMap(useForFitting,1))
        for j=1:length(expTime)
            if j==1
                fprintf(fid, '%s,%.1f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n',...
                              simData.DataNames{i},...
                              expTime(j),...
                              experimentalResults(j,i),...
                              experimentalResults(j,i),...
                              predictedResults(j,i),...
                              mseDRUG(i),...
                              maeDRUG(i),...
                              corrDRUG(i));
            else
                fprintf(fid, '%s,%.1f,%.3f,%.3f,%.3f,,,\n',...
                              simData.DataNames{i},...
                              expTime(j),...
                              experimentalResults(j,i),...
                              experimentalResults(j,i),...
                              predictedResults(j,i));
            end
        end
    else
         for j=1:length(expTime)
            if j==1
                fprintf(fid, '%s,%.1f,%.3f,,%.3f,%.3f,%.3f,%.3f\n',...
                              simData.DataNames{i},...
                              expTime(j),...
                              experimentalResults(j,i),...
                              predictedResults(j,i),...
                              mseDRUG(i),...
                              maeDRUG(i),...
                              corrDRUG(i));
            else
                fprintf(fid, '%s,%.1f,%.3f,,%.3f,,,\n',...
                              simData.DataNames{i},...
                              expTime(j),...
                              experimentalResults(j,i),...
                              predictedResults(j,i));
            end
         end
    end
end

% calculate bacterial contribution to serum metabolite
set(configsetObj, 'StopTime', max(expTime)+1);
simData = sbiosimulate(modelGutUniversal);
predictedResults = simData.Data;

if ismember('M_serum', simData.dataNames) && ismember('M_serumbact', simData.dataNames)
    metProfilesTime = simData.Time;
    metProfiles1 = predictedResults(:,ismember(simData.dataNames, 'M_serum'))/...
                    curVolumes(ismember(simData.dataNames, 'M_serum'));
    metProfiles2 = predictedResults(:,ismember(simData.dataNames, 'M_serumbact'))/...
                    curVolumes(ismember(simData.dataNames, 'M_serum'));

    hostContr = trapz(metProfilesTime,...
                      metProfiles1);
    bactContr = trapz(metProfilesTime,...
                      metProfiles2);
    bactFraction = bactContr/(bactContr+hostContr);

    fprintf(fid, 'Host contribution to M_serum = %.2f\n', hostContr);
    fprintf(fid, 'Bacterial contribution to M_serum = %.2f\n', bactContr);
    fprintf(fid, 'Bacterial fraction in contribution to M_serum = %.2f\n', bactFraction);
end
if ismember('M1_serum', simData.dataNames) && ismember('M1_serumbact', simData.dataNames)
    metProfilesTime = simData.Time;
    metProfiles1 = predictedResults(:,ismember(simData.dataNames, 'M1_serum'))/...
                    curVolumes(ismember(simData.dataNames, 'M1_serum'));
    metProfiles2 = predictedResults(:,ismember(simData.dataNames, 'M1_serumbact'))/...
                    curVolumes(ismember(simData.dataNames, 'M1_serum'));

    hostContr = trapz(metProfilesTime,...
                      metProfiles1);
    bactContr = trapz(metProfilesTime,...
                      metProfiles2);
    bactFraction = bactContr/(bactContr+hostContr);

    fprintf(fid, 'Host contribution to M1_serum = %.2f\n', hostContr);
    fprintf(fid, 'Bacterial contribution to M1_serum = %.2f\n', bactContr);
    fprintf(fid, 'Bacterial fraction in contribution to M1_serum = %.2f\n', bactFraction);
end
if ismember('M2_serum', simData.dataNames) && ismember('M2_serumbact', simData.dataNames)
    metProfilesTime = simData.Time;
    metProfiles1 = predictedResults(:,ismember(simData.dataNames, 'M2_serum'))/...
                    curVolumes(ismember(simData.dataNames, 'M2_serum'));
    metProfiles2 = predictedResults(:,ismember(simData.dataNames, 'M2_serumbact'))/...
                    curVolumes(ismember(simData.dataNames, 'M2_serum'));

    hostContr = trapz(metProfilesTime,...
                      metProfiles1);
    bactContr = trapz(metProfilesTime,...
                      metProfiles2);
    bactFraction = bactContr/(bactContr+hostContr);

    fprintf(fid, 'Host contribution to M2_serum = %.2f\n', hostContr);
    fprintf(fid, 'Bacterial contribution to M2_serum = %.2f\n', bactContr);
    fprintf(fid, 'Bacterial fraction in contribution to M2_serum = %.2f\n', bactFraction);
end
if ismember('M3_serum', simData.dataNames) && ismember('M3_serumbact',simData.dataNames)
    metProfilesTime = simData.Time;
    metProfiles1 = predictedResults(:,ismember(simData.dataNames, 'M3_serum'))/...
                    curVolumes(ismember(simData.dataNames, 'M3_serum'));
    metProfiles2 = predictedResults(:,ismember(simData.dataNames, 'M3_serumbact'))/...
                    curVolumes(ismember(simData.dataNames, 'M3_serum'));

    hostContr = trapz(metProfilesTime,...
                      metProfiles1);
    bactContr = trapz(metProfilesTime,...
                      metProfiles2);
    bactFraction = bactContr/(bactContr+hostContr);

    fprintf(fid, 'Host contribution to M3_serum = %.2f\n', hostContr);
    fprintf(fid, 'Bacterial contribution to M3_serum = %.2f\n', bactContr);
    fprintf(fid, 'Bacterial fraction in contribution to M3_serum = %.2f\n', bactFraction);
end
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
