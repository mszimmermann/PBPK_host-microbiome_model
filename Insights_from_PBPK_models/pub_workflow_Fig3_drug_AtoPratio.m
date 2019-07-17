%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script performs analysis of model parameters to study microbiome 
% effect on the systemic drug metabolite exposure for highly absorbed drugs.
% A. Fraction of bacterial contribution calculated for parameter sets that
% differ only by the percent of drug absorbed into the central compartment, 
% simulated as the ratio between kaP and kp1. Each point represents bacterial
% contribution to systemic drug metabolite exposure simulated for two 
% different absorption to propagation ratios 
% (kaP/kp1 = 0.1 and kaP/kp1 = 10, which corresponds to ~21% and ~93% of 
% initial drug dose absorbed into circulation, on x-axis and on y-axis)
% for a given parameter set. 
% B. Local sensitivity analysis of the fully parameterized pharmacokinetic 
% model for bacterial contribution to systemic drug metabolite exposure. 
% Drug absorption was kept high for all simulations (kaP to kp1 ratio = 10,
% which corresponds to 93% of the initial drug dose absorbed into circulation).
% Each point represents the calculated exposure to the drug for a given
% simulation, and color corresponds to the value of the tested parameter 
% in this simulation. 
% C. Bacterial contribution to systemic drug metabolite exposure under 
% different drug absorption to propagation ratios and bacterial to host 
% drug metabolism ratios. All the other parameters were kept the same for 
% each comparison. Biliary secretion parameter kEH was set to zero for all 
% simulations in C.  
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
% run the utility script to estimate how much of the initial drug dose
% was absorbed into systemic circulation at different 
% small-intestinal absorption to intestinal propagation ratios
AtoPratio_percentDrug = pub_get_AtoP_perent_absorption();
% define ratio to test for the second part (B)
testRatio = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the same model with different Absorption to Propagation ratios
% and test exposure to drug metabolite and bacterial contribution
% define two values of A to P ratios to test with each parameter set
Frange = [0.1 10];
% get indices of parameters of interest
% drug absorption from distal small intetsine
Fparam0 = ismember(baseModelParameterNames, 'k_aSI_P');
% initial drug absorption from small intestine
Fparam = ismember(baseModelParameterNames, 'k_a_P');
% drug propagation through small intestine
Pparam = ismember(baseModelParameterNames, 'k_p1');
metProfilesSimDataF = cell(nsim,length(Frange));
for f_i = 1:length(Frange)
    tic
    for sim_i = 1:nsim
        % set each parameter of the model to the current random value
        for pari = 1:length(baseModelParameterNames)
            newModel.Parameters(baseModelParameterIDX(pari)).Value = ...
                randParValues(pari, sim_i);
        end
        % set the ratio of A to P to the tested value f_i
         newModel.Parameters(baseModelParameterIDX(Fparam)).Value = ...
                Frange(f_i)*newModel.Parameters(baseModelParameterIDX(Pparam)).Value ;
         %% Simulate
         configsetObj = getconfigset(newModel);
         configsetObj.SolverOptions.SensitivityAnalysis = false;
         set(configsetObj, 'StopTime', 20);
         % get simulated data
         simData = sbiosimulate(newModel);
         % save the simulated data into a cell array
         metProfilesSimDataF{sim_i,f_i} = simData;
         
    end
    toctime = toc;
    fprintf('Finished iteration %d in %.3f sec\n', f_i, toctime)
end
% define names of metabolites of interest
% to calculate exposure for each simulated data
metProfilesMetF = cell(size(Frange));
metProfilesNames = {'M_serum', 'M_serumbact'};
n_metabolites = length(metProfilesNames);
for f_i = 1:length(Frange)
    metProfilesMet = zeros(nsim,n_metabolites);
    for sim_i = 1:nsim
        simData = metProfilesSimDataF{sim_i,f_i};
        for met_i = 1:n_metabolites
            metProfilesMet(sim_i,met_i) = trapz(simData.Time,...
                                               simData.Data(:,ismember(simData.dataNames,...
                                                                metProfilesNames{met_i})));
        end
    end
    metProfilesMetF{f_i} = metProfilesMet;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = figure('units','normalized','outerposition',[0 0 1 1]);
set(groot, 'defaultTextInterpreter','none')
set(groot, 'defaultAxesTickLabelInterpreter','none')
% combine M_host and M_bact from serum to assess total metabolite exposure
% in the serum
subplot(2,2,1);
scatter(metProfilesMetF{1}(:,2)./(metProfilesMetF{1}(:,2)+metProfilesMetF{1}(:,1)),...
        metProfilesMetF{end}(:,2)./(metProfilesMetF{end}(:,2)+metProfilesMetF{end}(:,1)),...
        15,'ko','filled');
hold on 
plot([0 1], [0 1], '--k')
axis square
xlabel({'Bacterial fraction,',sprintf('A/P = %.2f (%d%% drug absorbed)', Frange(1),...
                round(100*AtoPratio_percentDrug(AtoPratio_percentDrug(:,1)==Frange(1),2)))});
ylabel({'Bacterial fraction,',sprintf('A/P = %.2d (%d%% drug absorbed)', Frange(end),...
                round(100*AtoPratio_percentDrug(AtoPratio_percentDrug(:,1)==Frange(end),2)))});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for each parameter, fix it to the original value and alter the others
% to test which parameters influence bacterial contribution to systemic
% drug metabolite exposure under the condition of high drug absorption
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the model with sampled parameters locally
% for each parameter, fix it to the original value and alter the others
for par_i=1:length(baseModelParameterNames)
        newModel.Parameters(baseModelParameterIDX(par_i)).Value = ...
            baseModel.Parameters(baseModelParameterIDX(par_i)).Value;
end

metProfilesFtoPSimData = cell(length(baseModelParameterNames),nsim);
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
         % set the A to P ratio to the testRatio defined in the beginning 
         newModel.Parameters(baseModelParameterIDX(Fparam)).Value = ...
                testRatio*newModel.Parameters(baseModelParameterIDX(Pparam)).Value ;
         %% Simulate
         configsetObj = getconfigset(newModel);
         configsetObj.SolverOptions.SensitivityAnalysis = false;
         set(configsetObj, 'StopTime', 20);
         % get simulated data
         simData = sbiosimulate(newModel);
         % save simulated data to a cell array
         metProfilesFtoPSimData{par_i,sim_i} = simData;
    end
    toctime = toc;
    fprintf('Finished parameter %d in %.3f sec\n', par_i, toctime)    
end
% calculate exposure to the metabolites of interest (metProfilesNames)
n_metabolites = length(metProfilesNames);
metProfilesMet = cell(size(metProfilesNames));
for met_i = 1:n_metabolites
     curprofiles = zeros(length(baseModelParameterNames),nsim);
     for par_i = 1:length(baseModelParameterNames)   
        for sim_i = 1:nsim
         simData = metProfilesFtoPSimData{par_i,sim_i};
         curprofiles(par_i,sim_i) = trapz(simData.Time,...
                                          simData.Data(:,ismember(simData.dataNames,...
                                                         metProfilesNames{met_i})));
        end
     end
    metProfilesMet{met_i} = curprofiles;
end

%plot tornado plot of parameter effect on bacterial contribution to the
%systemic drug metabolite exposure
ax2 = subplot(2,2,[2, 4]);
hold on
for par_i = 1:length(baseModelParameterNames)
    scatter( metProfilesMet{2}(par_i,:)./...
             (metProfilesMet{1}(par_i,:)+metProfilesMet{2}(par_i,:)),...
             par_i*ones(nsim,1),...
             dotsize, log10(randParValues(par_i,:)), 'filled')%,'k')
end
title(sprintf('Bacterial contribution at A/P=%d (%d%% drug absorbed)',...
              testRatio,...
              round(100*AtoPratio_percentDrug(AtoPratio_percentDrug(:,1)==testRatio,2))),...
              'interpreter', 'none')
set(gca, 'YTick', 1:length(baseModelParameterNames))
set(gca, 'YTickLabel', baseModelParameterNames)
ylim([0 length(baseModelParameterNames)+1])
xlabel('Bacterial fraction in M_serum')
colormap(ax2,parula)
h = colorbar;
ylabel(h, 'Parameter value, log10')
set(gca,'Ydir','reverse')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate different A to P ratios and Bact to Host conversion ratios 
% when enterohepatic cycling = 0
% define range for A to P ratios
Frange = [0.1 1:10 100];
% define range for bacterial to host drug to metabolite conversion coefficients
BHrange = [0.1 1:10 100];
% get indeces of parameters of interest:
% initial drug absorption from small intestine
Fparam = ismember(baseModelParameterNames, 'k_a_P');
% drug absorption from distal small intestine
SIparam = ismember(baseModelParameterNames, 'k_aSI_P');
% enterohepatic transition from central compartment to galbladder 
EHparam = ismember(baseModelParameterNames, 'k_eh');
% drug propagation through the small intestine
Pparam = ismember(baseModelParameterNames, 'k_p1');
% drug to metabolite conversion by the host
Hparam = ismember(baseModelParameterNames, 'k_cH_M');
% drug to metabolite conversion by bacteria
Bparam = ismember(baseModelParameterNames, 'k_cB_M');

for pari = 1:length(baseModelParameterNames)
    newModel.Parameters(baseModelParameterIDX(pari)).Value = ...
        baseModel.Parameters(baseModelParameterIDX(pari)).Value;
end
newModel.Parameters(baseModelParameterIDX(EHparam)).Value = 0;
newModel.Parameters(baseModelParameterIDX(SIparam)).Value = 0;

metProfilesSimDataF = cell(length(Frange),length(BHrange));
for f_i = 1:length(Frange)
    for eh_i = 1:length(BHrange)
        tic
         % set the ratio of A to P to the current value f_i
         newModel.Parameters(baseModelParameterIDX(Fparam)).Value = ...
                Frange(f_i)*newModel.Parameters(baseModelParameterIDX(Pparam)).Value ;
         % set the ratio of bacterial to host drug conversion to the
         % current value eh_i
         newModel.Parameters(baseModelParameterIDX(Bparam)).Value = ...
                BHrange(eh_i)*newModel.Parameters(baseModelParameterIDX(Hparam)).Value ;
         
          %% Simulate
         configsetObj = getconfigset(newModel);
         configsetObj.SolverOptions.SensitivityAnalysis = false;
         set(configsetObj, 'StopTime', 20);
         % get simulated data
         simData = sbiosimulate(newModel);
         % save simulated data to a cell array
         metProfilesSimDataF{f_i, eh_i} = simData;
         toctime = toc;
         fprintf('Finished iteration %d %d in %.3f sec\n', f_i, eh_i, toctime)
    end
end
% define metabolites of interest to calculate exposure
metProfilesNames = {'M_serum', 'M_serumbact'};
n_metabolites = length(metProfilesNames);
metProfilesMet = cell(size(metProfilesNames));
for met_i = 1:n_metabolites
     curprofiles = zeros(size(metProfilesSimDataF,1),size(metProfilesSimDataF,2));
     for i = 1:size(metProfilesSimDataF,1)   
        for j = 1:size(metProfilesSimDataF,2)
            simData = metProfilesSimDataF{i,j};
            curprofiles(i,j) = trapz(simData.Time,...
                                     simData.Data(:,ismember(simData.dataNames,...
                                                         metProfilesNames{met_i})));
        end
     end
    metProfilesMet{met_i} = curprofiles;
end
% calculate bacterial contribution to systemic drug metabolite
metProfilesMetRatio = metProfilesMet{2}./...
                     (metProfilesMet{1}+metProfilesMet{2});
% convert A to P ratio to the percent of drug absorbed calculated in the beginning             
Frange_percent = cell(size(Frange));
for i=1:length(Frange)
    Frange_percent{i} = strcat(num2str(Frange(i)),...
                               sprintf(' (%d%% absorbed)',...
                                       round(100*AtoPratio_percentDrug(AtoPratio_percentDrug(:,1)==Frange(i),2))));
end

% plot heatmap of bacterial contribution at different combinations of 
% A to P ratio and bacterial to host drug convertion ratio
ax3 = subplot(2,2,3);
imagesc(metProfilesMetRatio)
set(gca, 'XTick', 1:length(BHrange))
set(gca, 'XTickLabel', BHrange)
set(gca, 'YTick', 1:length(Frange))
set(gca, 'YTickLabel', Frange_percent)
set(gca,'Ydir','normal')
xlabel('cB to cH ratio')
ylabel('Absorption to propagation ratio')
h = colorbar;
colormap(ax3, flipud(gray))
ylabel(h, 'Bacterial contribution to M_serum')
axis square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%