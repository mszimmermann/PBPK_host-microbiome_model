%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs analysis of model parameters to study the interplay
% between enterohepatic circulation and microbiome contribution to drug metabolism.
% B. Local sensitivity analysis of the fully parameterized pharmacokinetic 
% model for the difference of drug and drug metabolite exposure in serum 
% between two conditions with low and high biliary secretion (kEH). 
% C. Local sensitivity analysis of the fully parameterized pharmacokinetic 
% model for the difference of the drug and drug metabolite exposure in serum
% between two conditions with low and high bacterial glucuronidase 
% activity (kdglB). Each point represents the calculated exposure to
% the drug or metabolite for a given simulation, and color corresponds to 
% the value of the tested parameter in this simulation. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the parent folder with utility scripts and data folder
parentFolder =  strsplit(pwd, filesep);
parentFolder = strjoin(parentFolder(1:end-1), filesep);
addpath(genpath(parentFolder));
dataFolder = ['Data' filesep];
% define filename containing the model description 
modelfilename = [dataFolder 'example_extended_model_HYP_GLCRN.csv'];
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define names of parameters to change:
% k_eh - enterohepatic propagation from the cnetral compartment to
% galbladder and from galbladder to small intestine;
% k_dglB - bacterial coefficient of drug deglucuronidation
changeParamNames = {'k_eh', 'k_dglB'};
% define minimal and maximal values for the tested parameters
minParValue = 0;
maxParValue = 10;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each parameter, fix it to the original value and alter the others
fig = figure('units','normalized','outerposition',[0 0 1 1]);
spidx = 1; % index of subplots on the figure
for par_name = 1:length(changeParamNames)
    % fine the index of parameter to be tested
    curParamName = changeParamNames(par_name);
    changeParIdx = find(ismember(baseModelParameterNames,curParamName ));
    % save differences in exposure to drug and metabolite in serum and cecum
    metProfilesParentSerumDiff = zeros(length(baseModelParameterNames),nsim);
    metProfilesMetSerumDiff = zeros(length(baseModelParameterNames),nsim);
    metProfilesParentCecumDiff = zeros(length(baseModelParameterNames),nsim);
    metProfilesMetCecumDiff = zeros(length(baseModelParameterNames),nsim);
    % set the model parameter values to initial values
    for pari = 1:length(baseModelParameterNames)
            newModel.Parameters(baseModelParameterIDX(pari)).Value = ...
                baseModel.Parameters(baseModelParameterIDX(pari)).Value;
    end
    % get the serum volume of dustribution from the model
    for par_i = 1:length(baseModelParameterNames)
        tic;
        % reset the parameter value from the previous iteration to the initial value
        if par_i>1
            newModel.Parameters(baseModelParameterIDX(par_i-1)).Value = ...
                baseModel.Parameters(baseModelParameterIDX(par_i-1)).Value;
        end
        for sim_i = 1:nsim
            % set the value of the current tested parameter to random value
            newModel.Parameters(baseModelParameterIDX(par_i)).Value = ...
                    randParValues(par_i, sim_i);
            % set the value of the parameter of interest to min value
            newModel.Parameters(baseModelParameterIDX(changeParIdx)).Value = minParValue;
             %% Simulate
             configsetObj = getconfigset(newModel);
             configsetObj.SolverOptions.SensitivityAnalysis = false;
             set(configsetObj, 'StopTime', 20);
             % get simulated data
             simData = sbiosimulate(newModel);
             % calculate exposure to drug and metabolite in different tissues
             metProfilesParentSer0 = trapz(simData.Time,...
                                  simData.Data(:,ismember(simData.dataNames, 'P_serum')));
             metProfilesParentCec0 = trapz(simData.Time,...
                                  simData.Data(:,ismember(simData.dataNames, 'P_cecum')));
             metProfilesMetCec0 = trapz(simData.Time,...
                                  simData.Data(:,ismember(simData.dataNames, 'M_cecum')));
             metProfilesMetSer0 = trapz(simData.Time,...
                                  simData.Data(:,ismember(simData.dataNames, 'M_serum')))+...
                                  trapz(simData.Time,...
                                  simData.Data(:,ismember(simData.dataNames, 'M_serumbact')));
            % set the value of the parameter of interest to max value
             newModel.Parameters(baseModelParameterIDX(changeParIdx)).Value = maxParValue;
             %% Simulate
             configsetObj = getconfigset(newModel);
             configsetObj.SolverOptions.SensitivityAnalysis = false;
             set(configsetObj, 'StopTime', 20);
             % get simulated data
             simData = sbiosimulate(newModel);
            % calculate exposure to drug and metabolite in different tissues 
             metProfilesParentSer10 = trapz(simData.Time,...
                                  simData.Data(:,ismember(simData.dataNames, 'P_serum')));
             metProfilesParentCec10 = trapz(simData.Time,...
                                  simData.Data(:,ismember(simData.dataNames, 'P_cecum')));
             metProfilesMetCec10 = trapz(simData.Time,...
                                  simData.Data(:,ismember(simData.dataNames, 'M_cecum')));
             metProfilesMetSer10 = trapz(simData.Time,...
                                  simData.Data(:,ismember(simData.dataNames, 'M_serum')))+...
                                  trapz(simData.Time,...
                                  simData.Data(:,ismember(simData.dataNames, 'M_serumbact')));
            % calculate the difference of exposure to drug and metabolite
            % between the two paramater values
             metProfilesParentSerumDiff(par_i,sim_i) = metProfilesParentSer0-metProfilesParentSer10;
             metProfilesParentCecumDiff(par_i,sim_i) = metProfilesParentCec0-metProfilesParentCec10;
             metProfilesMetSerumDiff(par_i,sim_i) = metProfilesMetSer0-metProfilesMetSer10;
             metProfilesMetCecumDiff(par_i,sim_i) = metProfilesMetCec0-metProfilesMetCec10;
        end
        toctime = toc;
        fprintf('Finished parameter %d in %.3f sec\n', par_i, toctime)    
    end

    % plot tornado plots depicting effect of different parameters on
    % exposure difference to drug and drug metabolite
    for i=1:2
        subplot(1,4,spidx)
        hold on
        spidx = spidx+1;

        set(groot, 'defaultAxesTickLabelInterpreter','none')
        set(groot, 'defaultTextInterpreter','none')
        hold on
        switch i
            case 1           
                for par_i = 1:length(baseModelParameterNames)
                    scatter(metProfilesParentSerumDiff(par_i, :), par_i*ones(nsim,1),...
                        dotsize, log10(randParValues(par_i,:)), 'filled')%,'k')
                end
                title(['Pserum difference between 0 and 10 ' curParamName], 'interpreter', 'none')
            case 2
                for par_i = 1:length(baseModelParameterNames)
                    scatter(metProfilesMetSerumDiff(par_i, :), par_i*ones(nsim,1),...
                        dotsize, log10(randParValues(par_i,:)), 'filled')%,'k')
                end
                title(['Mserum difference between 0 and 10 ' curParamName], 'interpreter', 'none')
            case 3
                for par_i = 1:length(baseModelParameterNames)
                    scatter(metProfilesParentCecumDiff(par_i, :), par_i*ones(nsim,1),...
                        dotsize, log10(randParValues(par_i,:)), 'filled')%,'k')
                end
                title(['Pcecum difference between 0 and 10 ' curParamName], 'interpreter', 'none')
            case 4
                for par_i = 1:length(baseModelParameterNames)
                    scatter(metProfilesMetCecumDiff(par_i, :), par_i*ones(nsim,1),...
                        dotsize, log10(randParValues(par_i,:)), 'filled')%,'k')
                end
                title(['Mcecum difference between 0 and 10 ' curParamName], 'interpreter', 'none')            
        end
        set(gca, 'YTick', 1:length(baseModelParameterNames))
        set(gca, 'YTickLabel', baseModelParameterNames)
        ylim([0 length(baseModelParameterNames)+1])
        xlabel('Difference in serum AUC')
        h = colorbar;
        ylabel(h, 'Parameter value, log10')
        set(gca,'Ydir','reverse')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


