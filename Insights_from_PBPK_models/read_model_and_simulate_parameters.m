function [baseModel, newModel,...
          baseModelParameterNames,...
         baseModelParameterIDX,...
         randParValues] = read_model_and_simulate_parameters(modelfilename,...
                                                             nsim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
change_bioavailabilityF_flag = 0;
change_dose_flag = 0;
simulate_within_estimated_flag = 0;
plot_parameter_distribution_flag = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read the model from file
[baseModel] = create_model_from_file(modelfilename);
newModel = copyobj(baseModel);

% sensitivity analysis: simulate model with changing parameters
% get model parameter names
baseModelParameterNames = cell(length(baseModel.Parameters),1);
baseModelParameterIDX = zeros(length(baseModel.Parameters),1);
modelParameterEstimetesMEAN = zeros(length(baseModel.Parameters),1);

Fidx = 0;
Didx = 0;
for i=1:length(baseModel.Parameters)
    % only add parameters that are not constants for random sampling
    if ~(ismember(baseModel.Parameters(i).Tag, {'constant'}))
        baseModelParameterNames{i} = baseModel.Parameters(i).Name;
        modelParameterEstimetesMEAN(i) = baseModel.Parameters(i).Value;
        baseModelParameterIDX(i) = i;
    end
    if isequal(baseModel.Parameters(i).Name, 'F')
        Fidx = i;
    end
    if isequal(baseModel.Parameters(i).Name, 'D')
        Didx = i;
    end
end
% add drug dose D to simulated parameters if specified in the flag
if change_dose_flag
    baseModelParameterNames = [baseModelParameterNames; 'D'];
    baseModelParameterIDX(end+1) = Didx;
    modelParameterEstimetesMEAN(end+1) =  baseModel.Parameters(Didx).Value;
end
% add bioavailability F to simulated parameters if specified in the flag
if change_bioavailabilityF_flag
    baseModelParameterNames = [baseModelParameterNames; 'F'];
    baseModelParameterIDX(end+1) = Fidx;
    modelParameterEstimetesMEAN(end+1) =  baseModel.Parameters(Fidx).Value;
end

baseModelParameterIDX(cellfun(@(x) isempty(x), baseModelParameterNames))=[];
modelParameterEstimetesMEAN(cellfun(@(x) isempty(x), baseModelParameterNames))=[];
baseModelParameterNames(cellfun(@(x) isempty(x), baseModelParameterNames))=[];

% sort parameter names in specific order
parameterOrder = {'k_a_P','k_aSI_P','k_e_P',...
                  'k_p1','k_p2','k_p3','k_p4',...
                  'k_eh','k_aLI1_P',...
                  'k_cH_M','k_e_M',...
                  'k_cB_M', 'k_aLI1_M', ...
                  'k_gl_P','k_e_glP','k_dglB'};
[~, ~, idx] = intersect(parameterOrder, baseModelParameterNames, 'stable');
[~,idx_rest] = setdiff( baseModelParameterNames, baseModelParameterNames(idx));
idx = [idx; idx_rest];
baseModelParameterNames = baseModelParameterNames(idx);
baseModelParameterIDX = baseModelParameterIDX(idx);
modelParameterEstimetesMEAN = modelParameterEstimetesMEAN(idx);

% get the means and STD of the parameter estimates 
modelParameterEstimetesSTD = modelParameterEstimetesMEAN/2;

% get the intestinal propagation parameters k_p1-k_p5
propagationParameters = find(cellfun(@(x) ismember(x,...
                                 arrayfun(@(x) strcat('k_p',num2str(x)),[1:5], 'unif', 0)),...
                                 baseModelParameterNames));
% get the bioavailability parameter F
bioavailParameters = find(ismember(baseModelParameterNames, 'F'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if simulate within estimated range, use normal distribution with 
% estimated mean and STD for sampling
randParValues = randn(length(baseModelParameterNames),nsim);

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
    randParValuesNorm01 = normrnd(0, 1, [length(baseModelParameterNames),nsim]);
    nfactor = 1;
    randParValuesNorm = randn(length(baseModelParameterNames),nsim);
    for i=1:length(baseModelParameterNames)
        randParValuesNorm(i,:) = abs(randParValues(i,:)*modelParameterEstimetesSTD(i)*nfactor +...
                             modelParameterEstimetesMEAN(i));
    end
    %randParValues = randParValuesUnif;
    randParValues = 10.^randParValuesNorm01;
    % set norm distribution for propagation parameters
    randParValues(propagationParameters,:) = randParValuesNorm(propagationParameters,:);
    % set U[0 1] distribution for bioavailability F
    randParValues(bioavailParameters,:) = rand(length(bioavailParameters),nsim);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot parameter distributions
if plot_parameter_distribution_flag
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    spnum = ceil(size(randParValues,1)/3);
    for i=1:size(randParValues,1)
        subplot(3,spnum,i)
        if simulate_within_estimated_flag
            hist((randParValues(i,:)))
            hold on
            errorbar(modelParameterEstimetesMEAN(i), 300, modelParameterEstimetesSTD(i), 'horizontal', 'LineWidth', 3)
        else
            if ismember(i, propagationParameters) ||...
               ismember(i, bioavailParameters)
                hist((randParValues(i,:)))
                hold on
                errorbar(modelParameterEstimetesMEAN(i), 300, modelParameterEstimetesSTD(i), 'horizontal', 'LineWidth', 3)
            else
                hist(log10(randParValues(i,:)))
                hold on
                bar(log10(modelParameterEstimetesMEAN(i)), 300)
            end
        end
        title(baseModelParameterNames{i})
    end

    if simulate_within_estimated_flag
        suptitle({sprintf('Random parameter sampling for %s, n=%d', baseModel.Name, nsim),...
              'Distribution: Normal with mean and std of the estimated values'})
        % print figure to file
        orient landscape
    else
        suptitle({sprintf('Random parameter sampling for %s, n=%d', baseModel.Name, nsim),...
              ['Distribution: log10 uniform for non-propagation parameters,',...
              ' Normal with mean=estimatedValue and 3*std for propagation parameters']})
          % print figure to file
        orient landscape
    end
end
