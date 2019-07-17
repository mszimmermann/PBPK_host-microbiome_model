%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs analysis of model parameters to study microbiome 
% effect on the drug and metabolite levels in large intestine and serum. 
% A. Local sensitivity analysis of the fully parameterized pharmacokinetic 
% model for total drug levels in serum and cecum. 
% Each point represents the calculated exposure to the drug for a given 
% simulation on the x-axis, and color corresponds to the value of the 
% tested parameter in this simulation. 
% B. Simulated parent drug profiles in cecum and serum for two parameter 
% sets at initial values that differ only by the bacterial drug metabolism
% coefficient kcB. 
% C. Local sensitivity analysis of the fully parameterized pharmacokinetic 
% model for the difference of the drug exposure between two conditions with
% low and high bacterial drug metabolism. 
% D. Local sensitivity analysis of the fully parameterized pharmacokinetic 
% model for total metabolite levels in serum and cecum. Each point 
% represents the calculated exposure to the drug metabolite for a given 
% simulation on the x axis, and color corresponds to the value of the 
% tested parameter in this simulation. 
% E. Simulated drug metabolite profiles in cecum and serum for two 
% parameter sets that differ only by the bacterial drug metabolism 
% coefficient kcB. The parameter sets are the same as in B. 
% F. Local sensitivity analysis of the fully parameterized pharmacokinetic
% model for the difference of the drug metabolite exposure between two 
% conditions with low and high bacterial drug metabolism.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the parent folder with utility scripts and data folder
parentFolder =  strsplit(pwd, filesep);
parentFolder = strjoin(parentFolder(1:end-1), filesep);
addpath(genpath(parentFolder));
dataFolder = ['Data' filesep];
% define filename containing the model description 
modelfilename = [dataFolder 'example_extended_model_HYP_GF_CV.csv'];
% define number of simulations
nsim = 100;
% read model from the file and simulate values for model parameters
[baseModel, newModel,...
         baseModelParameterNames,...
         baseModelParameterIDX,...
         randParValues] = read_model_and_simulate_parameters(modelfilename,...
                                                             nsim);
% define size of dot for plotting tornado plots
dotsize = 100;
% define the parameter to analyze 
% (in this script, bacterial contribution to systemic drug and metabolite
% levels will be tested with varying k_cB_M parameter of bacterial drug
% metabolism in the cecum)
par_name = 'k_cB_M';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the model with sampled parameters locally:
% for each parameter, fix it to the original value and alter the values of 
% all other parameters
% save simulated data for each parameter and simulation
metProfilesParameterSimData = cell(length(baseModelParameterNames),nsim);
for par_i = 1:length(baseModelParameterNames)
    tic;
    % reset the parameter value from the previous iteration to the initial value
    if par_i>1
        newModel.Parameters(baseModelParameterIDX(par_i-1)).Value = ...
            baseModel.Parameters(baseModelParameterIDX(par_i-1)).Value;
    end
     
    for sim_i = 1:nsim
        % set the value of the currently tested parameter to the current
        % random value
        newModel.Parameters(baseModelParameterIDX(par_i)).Value = ...
                randParValues(par_i, sim_i);
    
         %% Simulate
         configsetObj = getconfigset(newModel);
         configsetObj.SolverOptions.SensitivityAnalysis = false;
         set(configsetObj, 'StopTime', 20);
         % get simulated data
         simData = sbiosimulate(newModel);

         % save the simulated data into a cell array
         metProfilesParameterSimData{par_i,sim_i} = simData;
    end
    toctime = toc;
    fprintf('Finished parameter %d in %.3f sec\n', par_i, toctime)    
end
% define names of metabolites of interest
% to calculate exposure for each simulated data
metProfilesNames = {'M_serum', 'M_serumbact', 'P_serum', 'P_cecum', 'M_cecum'};
n_metabolites = length(metProfilesNames);
metProfilesMet = cell(size(metProfilesNames));
for met_i = 1:n_metabolites
     curprofiles = zeros(length(baseModelParameterNames),nsim);
     for par_i = 1:length(baseModelParameterNames)   
        for sim_i = 1:nsim
         simData = metProfilesParameterSimData{par_i,sim_i};
         curprofiles(par_i,sim_i) = trapz(simData.Time,...
                                          simData.Data(:,ismember(simData.dataNames,...
                                                         metProfilesNames{met_i})));
        end
     end
    metProfilesMet{met_i} = curprofiles;
end
% combine M_host and M_bact from serum to assess total metabolite exposure
% in the serum
metProfilesMet{end+1} = metProfilesMet{1}+metProfilesMet{2};
metProfilesNames{end+1} = 'M_serum_total';
% plot tornado plots for different parameters and metabolite exposures
fig = figure('units','normalized','outerposition',[0 0 1 1]);
spidx = 1;
for met_i = 3:length(metProfilesNames)
    subplot(2,2,spidx)
    spidx = spidx+1;
    set(groot, 'defaultTextInterpreter', 'none')
    set(groot, 'defaultAxesTickLabelInterpreter','none')
    hold on
    for par_i = 1:length(baseModelParameterNames)
        scatter(metProfilesMet{met_i}(par_i, :), par_i*ones(nsim,1),...
            dotsize, log10(randParValues(par_i,:)), 'filled')%,'k')
    end
    title(metProfilesNames{met_i}, 'interpreter', 'none')
    set(gca, 'YTick', 1:length(baseModelParameterNames))
    set(gca, 'YTickLabel', baseModelParameterNames)
    ylim([0 length(baseModelParameterNames)+1])
    xlabel('Exposure (AUC)')
    h = colorbar;
    ylabel(h, 'Parameter value, log10')
    set(gca,'Ydir','reverse')
    axis square
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot profiles of drug and metabolites in cecum and serum
% define colors for low and high bacterial metabolism parameters
mycolors = [0 115 179;...
            230 153 0]./255;
% select the parameter defined in the beginning of the script (par_name)
par_i = ismember(baseModelParameterNames, par_name);

% get the minimum and maximum values for comparison of extremes
[minval, sim1] = min(randParValues(par_i,:));
[maxval, sim2] = max(randParValues(par_i,:));

% plot figure for drug and metabolite exposure in cecum and serum at 
% min and max values of the tested parameter
fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(groot, 'defaultAxesTickLabelInterpreter','none')
set(groot, 'defaultTextInterpreter','none')
spidx = 1;
for met_i = 3:5
    subplot(2,2,spidx)
    spidx = spidx+1;
    hold on
    simData = metProfilesParameterSimData{par_i, sim1};
    plot(simData.Time,...
         simData.Data(:,ismember(simData.dataNames,metProfilesNames{met_i})), ...
         '-.', 'LineWidth', 2, 'Color', mycolors(1,:));
    simData = metProfilesParameterSimData{par_i, sim2};
    plot(simData.Time,...
         simData.Data(:,ismember(simData.dataNames,metProfilesNames{met_i})), ...
         '--', 'LineWidth', 2, 'Color', mycolors(2,:));
    legend({sprintf('%s=%.3f', par_name, minval),...
            sprintf('%s=%.3f', par_name, maxval)})
    title(metProfilesNames{met_i})
    xlabel('Time, h')
end

subplot(2,2,spidx)
met_i1 = 1;
met_i2 = 2;
hold on
simData = metProfilesParameterSimData{par_i, sim1};
plot(simData.Time,...
     simData.Data(:,ismember(simData.dataNames,metProfilesNames{met_i1}))+...
     simData.Data(:,ismember(simData.dataNames,metProfilesNames{met_i2})), ...
     '-.', 'LineWidth', 2, 'Color', mycolors(1,:));
exp1 = trapz(simData.Time,...
     simData.Data(:,ismember(simData.dataNames,metProfilesNames{met_i1}))+...
     simData.Data(:,ismember(simData.dataNames,metProfilesNames{met_i2})));
simData = metProfilesParameterSimData{par_i, sim2};
plot(simData.Time,...
     simData.Data(:,ismember(simData.dataNames,metProfilesNames{met_i1}))+...
     simData.Data(:,ismember(simData.dataNames,metProfilesNames{met_i2})), ...
     '--', 'LineWidth', 2, 'Color', mycolors(2,:));
exp2 = trapz(simData.Time,...
     simData.Data(:,ismember(simData.dataNames,metProfilesNames{met_i1}))+...
     simData.Data(:,ismember(simData.dataNames,metProfilesNames{met_i2}))); 

legend({sprintf('%s=%.3f', par_name, minval),...
            sprintf('%s=%.3f', par_name, maxval)})
title('M_serum total')
xlabel('Time, h')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% analyze which parameter affects most the difference in exposure to
% systemic drug and metabolite between the metabolite metabolism conditions
% For each parameter, fix it to the original value and alter the others
changeParIdx = find(ismember(baseModelParameterNames, par_name));
% save the difference between drug and metabolite exposure at differet
% parameter values (min and max)
metProfilesParentSerumDiff_M = zeros(length(baseModelParameterNames),nsim);
metProfilesParentSerumDiff_P = zeros(length(baseModelParameterNames),nsim);

for par_i = 1:length(baseModelParameterNames)
    tic;
  
    % reset the parameter value from the previous iteration to the initial value
    if par_i>1
        newModel.Parameters(baseModelParameterIDX(par_i-1)).Value = ...
            baseModel.Parameters(baseModelParameterIDX(par_i-1)).Value;
    end
        
    for sim_i = 1:nsim
        % set the value of the currently tested parameter to the current
        % random value
        newModel.Parameters(baseModelParameterIDX(par_i)).Value = ...
                randParValues(par_i, sim_i);
        % set the value of the parameter of interest to minimum value
        newModel.Parameters(baseModelParameterIDX(changeParIdx)).Value = minval;
         %% Simulate
         configsetObj = getconfigset(newModel);
         configsetObj.SolverOptions.SensitivityAnalysis = false;
         set(configsetObj, 'StopTime', 20);

         simData = sbiosimulate(newModel);
         % save drug and metabolite exposure values at min parameter value
         metProfilesP0 = trapz(simData.Time,...
                              simData.Data(:,ismember(simData.dataNames, 'P_serum')));
         metProfilesM0 = trapz(simData.Time,...
                              simData.Data(:,ismember(simData.dataNames, 'M_serum')))+...
                       trapz(simData.Time,...
                              simData.Data(:,ismember(simData.dataNames, 'M_serumbact')));
       
                          
         % set the value of the parameter of interest to maximum value
         newModel.Parameters(baseModelParameterIDX(changeParIdx)).Value = maxval;
         %% Simulate
         configsetObj = getconfigset(newModel);
         configsetObj.SolverOptions.SensitivityAnalysis = false;
         set(configsetObj, 'StopTime', 20);

         simData = sbiosimulate(newModel);
         % save drug and metabolite exposure values at max parameter value
         metProfilesP100 = trapz(simData.Time,...
                              simData.Data(:,ismember(simData.dataNames, 'P_serum')));
         metProfilesM100 = trapz(simData.Time,...
                              simData.Data(:,ismember(simData.dataNames, 'M_serum')))+...
                       trapz(simData.Time,...
                              simData.Data(:,ismember(simData.dataNames, 'M_serumbact')));
         % save the difference between drug and metabolite exposures
         % between the min and max parameter value conditions
         metProfilesParentSerumDiff_M(par_i,sim_i) = (metProfilesM0-metProfilesM100);
         metProfilesParentSerumDiff_P(par_i,sim_i) = (metProfilesP0-metProfilesP100);
    end
    toctime = toc;
    fprintf('Finished parameter %d in %.3f sec\n', par_i, toctime)    
end

% plot tornado plot for the differences in drug and metabolite systemic
% exposures at min and max values of parameter of interest
fig = figure('units','normalized','outerposition',[0 0 1 1]);
spidx = 1;
for i=1:2
    subplot(2,1,spidx)
    spidx = spidx+1;
    set(groot, 'defaultAxesTickLabelInterpreter','none')
    set(groot, 'defaultTextInterpreter','none')
    hold on
    for par_i = 1:length(baseModelParameterNames)
        if i==1
            scatter(metProfilesParentSerumDiff_P(par_i, 1:nsim), par_i*ones(nsim,1),...
                dotsize, log10(randParValues(par_i,1:nsim)), 'filled')%,'k')
        end
        if i==2
            scatter(metProfilesParentSerumDiff_M(par_i, 1:nsim), par_i*ones(nsim,1),...
                dotsize, log10(randParValues(par_i,1:nsim)), 'filled')%,'k')
        end
    end
    set(gca, 'YTick', 1:length(baseModelParameterNames))
    set(gca, 'YTickLabel', baseModelParameterNames)
    ylim([0 length(baseModelParameterNames)+1])
    h = colorbar;
    ylabel(h, 'Parameter value, log10')
    set(gca,'Ydir','reverse')
    if i==1
        xlabel('Difference in P_serum AUC')
        title(sprintf('Pserum difference between %.3f and %.3f k_c_bact',...
                       minval, maxval), 'interpreter', 'none')
    end
    if i==2
        xlabel('Difference in M_serum AUC')
        title(sprintf('Mserum difference between %.3f and %.3f k_c_bact',...
                      minval, maxval), 'interpreter', 'none')
    end
    axis square
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
