function perform_global_sensitivity_analysis(baseModel,...
                                             resultsHost,...
                                             resultsHostBact,...
                                             nsim,...
                                             change_bioavailabilityF_flag,...
                                             change_dose_flag,...
                                             simulate_within_estimated_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function persorms global sensitivity analysis 
% by randomly sampling form the parameter space
% and assessing the changes in 
% 1) total drug metabolite in serum
% 2) host-derived drug metabolite in serum
% 3) bacteria-derived drug metabolite in serum
% 4) relative bacterial contribution to drug metabolite in serum
% Input:
%     baseModel - model for which to perform sensitivity analysis
%     resultsHost - estimated parameters for the host-only model (used as
%                   initial values for sampling around them)
%     resultsHostBact - estimated parameters for the host-only model
%                   (used as initial values for sampling around them)
%     nsim - number of random parameter simulations
%     change_bioavailabilityF_flag - flag indicating whether to sample F
%     change_dose_flag - flag indicating whether to sample drug dose D
%     simulate_within_estimated_flag - flag indicating whether to simulate
%           within estimated parameter intervals (normal distibution);
%           if ==0, sample from uniform distribution between -3 and 3 
%           in log scale for all parameters apart from intestinal propagation
%           (sampled with normal mean +/- 3*std of the estimates)
%           and F (sampled uniform between 0 and 1).
% Output: 
%     tornado plots of correlation between parameters and 1)-4) values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                             
                                             
% sensitivity analysis: simulate model with changing parameters
% get model parameter anames
baseModelParameterNames = cell(length(baseModel.Parameters),1);
for i=1:length(baseModel.Parameters)
    % only add parameters that are not constants for random sampling
    if ~(ismember(baseModel.Parameters(i).Tag, {'constant'}))
        baseModelParameterNames{i} = baseModel.Parameters(i).Name;
    end
end
% add drug dose D to simulated parameters if specified in the flag
if change_dose_flag
    baseModelParameterNames = [baseModelParameterNames; 'D'];
end
% add bioavailability F to simulated parameters if specified in the flag
if change_bioavailabilityF_flag
    baseModelParameterNames = [baseModelParameterNames; 'F'];
end
baseModelParameterNames(cellfun(@(x) isempty(x), baseModelParameterNames))=[];
% get the intestinal propagation parameters k_p1-k_p5
propagationParameters = find(cellfun(@(x) ismember(x,...
                                 arrayfun(@(x) strcat('k_p',num2str(x)),[1:5], 'unif', 0)),...
                                 baseModelParameterNames));
% get the bioavailability parameter F
bioavailParameters = find(ismember(baseModelParameterNames, 'F'));
% get the means and STD of the parameter estimates 
modelParameterEstimetesMEAN = zeros(length(baseModel.Parameters),1);
modelParameterEstimetesSTD = zeros(length(baseModel.Parameters),1);

% extract information from the imput parameter estimate structures
for i=1:length(baseModelParameterNames)
    idx = ismember(table2cell(resultsHost.ParameterEstimates(:,'Name')), baseModelParameterNames{i});
    if nnz(idx)
        modelParameterEstimetesMEAN(i) = resultsHost.ParameterEstimates{idx,'Estimate'};
        modelParameterEstimetesSTD(i) = resultsHost.ParameterEstimates{idx,'StandardError'};
    end
    idx = ismember(table2cell(resultsHostBact.ParameterEstimates(:,'Name')), baseModelParameterNames{i});
    if nnz(idx)
        modelParameterEstimetesMEAN(i) = resultsHostBact.ParameterEstimates{idx,'Estimate'};
        modelParameterEstimetesSTD(i) = resultsHostBact.ParameterEstimates{idx,'StandardError'};
    end
    % check whether the parameter is defined in the input file (Notes==0)
    if str2double(baseModel.Parameters(i).Notes) == 0
        modelParameterEstimetesMEAN(i) = baseModel.Parameters(i).Value;
        % heuristic STD - half of the mean
        modelParameterEstimetesSTD(i) = baseModel.Parameters(i).Value/2;
    end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
randParValues = randn(length(baseModelParameterNames),nsim);
% if simulate within estimated range, use normal distribution with 
% estimated mean and STD for sampling
if simulate_within_estimated_flag
    nfactor = 1;
    randParValuesNorm = randn(length(baseModelParameterNames),nsim);
    for i=1:length(baseModelParameterNames)
        randParValuesNorm(i,:) = abs(randParValues(i,:)*modelParameterEstimetesSTD(i)*nfactor +...
                             modelParameterEstimetesMEAN(i));
    end
    randParValues=randParValuesNorm;
else
% if simulate a wide range, use
% log-uniform U[-3 3] distribution for all non-propagation parameters, 
% normal distributon N(mean, 3*std) for propagation parameters, 
% and U[0 1] for bioavailability F
    parRange = 10.^([-3:0.1:3]);
    randParValues = randn(length(baseModelParameterNames),nsim);
    randParIDX = randi(length(parRange), length(baseModelParameterNames),nsim);
    randParValuesUnif = arrayfun(@(x) parRange(x), randParIDX);

    nfactor = 3;
    randParValuesNorm = randn(length(baseModelParameterNames),nsim);
    for i=1:length(baseModelParameterNames)
        randParValuesNorm(i,:) = abs(randParValues(i,:)*modelParameterEstimetesSTD(i)*nfactor +...
                             modelParameterEstimetesMEAN(i));
    end
    randParValues = randParValuesUnif;
    % set norm distribution for propagation parameters
    randParValues(propagationParameters,:) = randParValuesNorm(propagationParameters,:);
    % set U[0 1] distribution for bioavailability F
    randParValues(bioavailParameters,:) = rand(length(bioavailParameters),nsim);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot parameter distributions
fig = figure('units','normalized','outerposition',[0 0 1 1]);
spnum = ceil(size(randParValues,1)/3);
for i=1:size(randParValues,1)
    subplot(3,spnum,i)
    if ismember(i, propagationParameters) ||...
       ismember(i, bioavailParameters)
        hist((randParValues(i,:)))
    else
        hist(log10(randParValues(i,:)))
    end
    title(baseModelParameterNames{i})
end

suptitle({sprintf('Random parameter sampling for %s, n=%d', baseModel.Name, nsim),...
          ['Distribution: log10 uniform for non-propagation parameters,',...
          ' Normal with mean=estimatedValue and 3*std for propagation parameters']})
% print figure to file
orient landscape
print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
      ['hist_parameter_sampling_for_', baseModel.Name, '.pdf'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot subsample of pairwise parameter distributions
fig = figure('units','normalized','outerposition',[0 0 1 1]);
idx = 1:length(baseModelParameterNames);
nselectsim = 50;
randidx500 = randi(nsim,nselectsim,1);
for i=1:length(baseModelParameterNames)
    for j=i:length(baseModelParameterNames)
        subplot(length(baseModelParameterNames)+1,length(baseModelParameterNames)+1,...
               (i-1)*(length(baseModelParameterNames)+1)+j)
        if i==j
            if ~ismember(idx(i), [propagationParameters; bioavailParameters])
                hist(log10(randParValues(idx(i),randidx500)))
            else
                hist(randParValues(idx(i),randidx500))
            end
        else
            if ~ismember(idx(i), [propagationParameters; bioavailParameters])
                randVal_i = log10(randParValues(idx(i),randidx500));
            else
                randVal_i = randParValues(idx(i),randidx500);
            end
            if ~ismember(idx(j), [propagationParameters; bioavailParameters])
                randVal_j = log10(randParValues(idx(j),randidx500));
            else
                randVal_j = randParValues(idx(j),randidx500);
            end
            plot(randVal_i, randVal_j, '.')
        end
        if i==1
            title(baseModelParameterNames{idx(j)})
        end
        if j==i
            ylabel(baseModelParameterNames{idx(i)})
        end
        set(gca, 'fontSize', 6)
    end
end
suptitle(sprintf('Random parameter sampling for %s, n=%d out of %d',...
                 baseModel.Name, nselectsim, nsim))
% print figure to file
orient landscape
print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
      ['scatter_pairwise_parameter_sampling_for_', baseModel.Name, '.pdf'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate the model with sampled parameters and calculate drug metabolite
% exposure and bacterial contribution
simData = sbiosimulate(baseModel);
nbact_metabolites = sum(cellfun(@(x) contains(x, '_serumbact'), simData.dataNames));

metProfilesParHost = zeros(nsim,nbact_metabolites);
metProfilesParBact = zeros(nsim,nbact_metabolites);
metProfilesNames = cell(1,nbact_metabolites);
% get the serum volume of dustribution from the model
VbloodDRUG = sbioselect(baseModel, 'Name','Vserum');
VbloodDRUG = VbloodDRUG.Value;
for sim_i = 1:nsim
    for pari = 1:length(baseModelParameterNames)
        baseModel.Parameters(pari).Value = randParValues(pari, sim_i);
    end
     %% Simulate
     configsetObj = getconfigset(baseModel);
     configsetObj.SolverOptions.SensitivityAnalysis = false;
     set(configsetObj, 'StopTime', 20);
     simData = sbiosimulate(baseModel);

     if nnz(simData.Data(:,ismember(simData.dataNames, 'M_serumbact'))>0)
         metProfilesParHost(sim_i) = trapz(simData.Time,...
                                           simData.Data(:,ismember(simData.dataNames, 'M_serum'))/...
                                           VbloodDRUG);
         metProfilesParBact(sim_i) = trapz(simData.Time,...
                                           simData.Data(:,ismember(simData.dataNames, 'M_serumbact'))/...
                                           VbloodDRUG);
         metProfilesNames(1) = {'M_serumbact'};
     end
     if nnz(simData.Data(:,ismember(simData.dataNames, 'M1_serumbact'))>0)
         metProfilesParHost(sim_i,1) = trapz(simData.Time,...
                                           simData.Data(:,ismember(simData.dataNames, 'M1_serum'))/...
                                           VbloodDRUG);
         metProfilesParBact(sim_i,1) = trapz(simData.Time,...
                                           simData.Data(:,ismember(simData.dataNames, 'M1_serumbact'))/...
                                           VbloodDRUG);
         metProfilesNames(1) = {'M1_serumbact'};
     end
     if nnz(simData.Data(:,ismember(simData.dataNames, 'M2_serumbact'))>0)
         metProfilesParHost(sim_i,2) = trapz(simData.Time,...
                                           simData.Data(:,ismember(simData.dataNames, 'M2_serum'))/...
                                           VbloodDRUG);
         metProfilesParBact(sim_i,2) = trapz(simData.Time,...
                                           simData.Data(:,ismember(simData.dataNames, 'M2_serumbact'))/...
                                           VbloodDRUG);
         metProfilesNames(2) = {'M2_serumbact'};
     end
     if nnz(simData.Data(:,ismember(simData.dataNames, 'M3_serumbact'))>0)
         metProfilesParHost(sim_i,3) = trapz(simData.Time,...
                                           simData.Data(:,ismember(simData.dataNames, 'M3_serum'))/...
                                           VbloodDRUG);
         metProfilesParBact(sim_i,3) = trapz(simData.Time,...
                                           simData.Data(:,ismember(simData.dataNames, 'M3_serumbact'))/...
                                           VbloodDRUG);
         metProfilesNames(3) = {'M3_serumbact'};
     end
end

metProfilesParTotal = metProfilesParHost+metProfilesParBact;
metContributionBact = metProfilesParBact./metProfilesParTotal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TORNADO PLOT
metProfilesParTotal_all = metProfilesParTotal;
metProfilesParHost_all = metProfilesParHost;
metProfilesParBact_all = metProfilesParBact;
metContributionBact_all = metContributionBact;

for met_i = 1:size(metProfilesParBact_all,2)
    metProfilesParTotal = metProfilesParTotal_all(:,met_i);
    metProfilesParHost = metProfilesParHost_all(:,met_i);
    metProfilesParBact = metProfilesParBact_all(:,met_i);
    metContributionBact = metContributionBact_all(:,met_i);
    % skip if metabolite is not produced
    if nnz(metProfilesParTotal)==0
        continue
    end
    corrTotal = partialcorri(metProfilesParTotal, randParValues');
    corrHost = partialcorri(metProfilesParHost, randParValues');
    corrBact = partialcorri(metProfilesParBact, randParValues');
    corrContrBact = partialcorri(metContributionBact, randParValues');

    parNames = baseModelParameterNames;

    [~, idx] = sort(abs(corrTotal), 'ascend');
    nancorr = corrTotal(idx);
    idx = [idx(isnan(nancorr)) idx(~isnan(nancorr)) ];

    corrTotal = corrTotal(idx);
    corrHost = corrHost(idx);
    corrBact = corrBact(idx);
    corrContrBact = corrContrBact(idx);
    parNames = parNames(idx);

    fig = figure('units','normalized','outerposition',[0 0 1 1]);

    subplot(1,4,1)
    barh(corrTotal)
    set(gca, 'YTick', 1:length(parNames))
    set(gca, 'YTickLabel', parNames)
    title('Total metabolite')
    xlim([-1 1])

    subplot(1,4,2)
    barh(corrHost)
    set(gca, 'YTick', 1:length(parNames))
    set(gca, 'YTickLabel', parNames)
    title('Host metabolite')
    xlim([-1 1])

    subplot(1,4,3)
    barh(corrBact)
    set(gca, 'YTick', 1:length(parNames))
    set(gca, 'YTickLabel', parNames)
    title('Bacterial metabolite')
    xlim([-1 1])

    subplot(1,4,4)
    barh(corrContrBact)
    set(gca, 'YTick', 1:length(parNames))
    set(gca, 'YTickLabel', parNames)
    title('Bacterial contribution')
    xlim([-1 1])

    suptitle(sprintf('%s %s parcorri partical correlation coefficients',...
             metProfilesNames{met_i},...
             baseModel.Name))
    orient landscape
    print(fig, '-painters', '-dpdf', '-r600', '-bestfit',...
          ['tornado_parcorri_parameter_sampling_for_',...
          baseModel.Name, '_',...
          metProfilesNames{met_i}, '.pdf'])
end