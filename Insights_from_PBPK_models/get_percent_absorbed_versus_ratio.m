%

nsim = 1000;
change_bioavailabilityF_flag = 0;
change_dose_flag = 0;
simulate_within_estimated_flag = 0;

modelfilename = ['Data\example_extended_model_HYP_SI_SERUM.csv'];
[baseModel] = create_model_from_file(modelfilename);
newModel = copyobj(baseModel);
configsetObj = getconfigset(newModel);
configsetObj.SolverOptions.SensitivityAnalysis = false;
set(configsetObj, 'StopTime', 20);
simData = sbiosimulate(newModel);
fig = figure('units','normalized','outerposition',[0 0 1 1]);
for i=1:length(simData.dataNames)
    subplot(3,ceil(length(simData.dataNames)/3),i)
    plot(simData.Time,simData.Data(:,i), 'LineWidth', 2);
    title(simData.DataNames{i})
    xlim([0 20])
    axis square
end

Frange = [0 0.01 0.1 1:10 20:100];
Fparam = ismember(baseModelParameterNames, 'k_a_P');
Pparam = ismember(baseModelParameterNames, 'k_p1');
% Fdrug = baseF;
% dose = baseD;
% resetModelParameters;

for pari = 1:length(baseModelParameterNames)
    newModel.Parameters(baseModelParameterIDX(pari)).Value = ...
        baseModel.Parameters(baseModelParameterIDX(pari)).Value;
end
newModel.Parameters(baseModelParameterIDX(ismember(baseModelParameterNames, 'k_aSI_P'))).Value = 0;
    
metProfilesSimDataF = cell(length(Frange),1);
for f_i = 1:length(Frange)
        tic


         newModel.Parameters(baseModelParameterIDX(Fparam)).Value = ...
                Frange(f_i)*newModel.Parameters(baseModelParameterIDX(Pparam)).Value ;
         
          %% Simulate
         configsetObj = getconfigset(newModel);
         configsetObj.SolverOptions.SensitivityAnalysis = false;
         set(configsetObj, 'StopTime', 20);

         simData = sbiosimulate(newModel);

         metProfilesSimDataF{f_i} = simData;
         toctime = toc;
         fprintf('Finished iteration %d in %.3f sec\n', f_i,toctime)
end

metProfilesNames = {'P_serum', 'P_si1'};
n_metabolites = length(metProfilesNames);
metProfilesMet = cell(size(metProfilesNames));
for met_i = 1:n_metabolites
     curprofiles = zeros(size(metProfilesSimDataF,1),size(metProfilesSimDataF,2));
     for i = 1:size(metProfilesSimDataF,1)   
        for j = 1:size(metProfilesSimDataF,2)
            simData = metProfilesSimDataF{i,j};
%             curprofiles(i,j) = trapz(simData.Time,...
%                                      simData.Data(:,ismember(simData.dataNames,...
%                                                          metProfilesNames{met_i})));
            curprofiles(i,j) = max(simData.Data(:,ismember(simData.dataNames,...
                                                         metProfilesNames{met_i})));
        end
     end
    metProfilesMet{met_i} = curprofiles;
end


