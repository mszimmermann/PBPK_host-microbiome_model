%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define files containing the models and data for first step:
% estimating host parameters (datafilenames1)
% and second step: estimating bacterial parameters (datafilenames2)
dataFolder = ['Data' filesep];
modelfilenames = {'example_basic_model_BRV_mut_wt.csv' ...
                  'example_basic_model_BRV_GF_CV.csv' ...
                  'example_basic_model_SRV_GF_CV.csv' ...
                  'example_basic_model_CLZ_GF_CV.csv' ...
                  'example_extended_model_CLZ_GF_CV.csv' ...
                  };
datafilenames1 = {'example_data_BRV_MUT.csv' ... 
                  'example_data_BRV_GF.csv'...
                  'example_data_SRV_GF.csv'...
                  'example_data_CLZ_GF.csv'...
                  'example_data_CLZ_GF_extended.csv'...
                  };
datafilenames2 = {'example_data_BRV_WT.csv'...
                  'example_data_BRV_CV.csv'...
                  'example_data_SRV_CV.csv'...
                  'example_data_CLZ_CV.csv'...
                  'example_data_CLZ_CV_extended.csv'
                   };
% flag indicating whether to perform sensitivity analysis
% (change to 1 to perfrom sensitivity analysis)
perform_sensitivity_analysis_flag = 0;              
% build the models
for files_i = 1:length(modelfilenames)
    modelfilename = [dataFolder modelfilenames{files_i}];
    datafilename1 = [dataFolder datafilenames1{files_i}];
    datafilename2 = [dataFolder datafilenames2{files_i}];
    % define model output files with out_ preffix and input file name
    outfilename1 = [dataFolder 'out_' datafilenames1{files_i}];
    outfilename2 = [dataFolder 'out_' datafilenames2{files_i}];
    
    % build model according to the definition in the table 
    [modelGutUniversal] = create_model_from_file(modelfilename);
    % load experimental data
    [t,metNamesMap, gd, useForFitting] = load_data_from_file(datafilename1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract model parameter names and initial values
    parameterNames = cell(size(modelGutUniversal.Parameters));
    parameterInitialValues = zeros(size(modelGutUniversal.Parameters));
    parameterNames_fitdata = zeros(size(modelGutUniversal.Parameters));
    for i=1:length(modelGutUniversal.Parameters)
        parameterNames{i} = modelGutUniversal.Parameters(i).Name;
        parameterNames_fitdata(i) = str2double(modelGutUniversal.Parameters(i).Notes);
        parameterInitialValues(i) = modelGutUniversal.Parameters(i).Value;
        if parameterNames_fitdata(i)~=0
            modelGutUniversal.Parameters(i).Value = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract volumes of the compartments for the first step (GF) and
    % second step (CVR)
    % Note: GF and CVR are historical names, if both mouse models are 
    % monocolonized, for example, the volumes for GF and CVR can be defined
    % the same
    VbloodDRUG = modelGutUniversal.Parameters(ismember(parameterNames, 'Vserum')).Value;
    VcecGF = modelGutUniversal.Parameters(ismember(parameterNames, 'Vcecum')).Value;
    VcolGF = modelGutUniversal.Parameters(ismember(parameterNames, 'Vcolon')).Value;
    VfecGF = modelGutUniversal.Parameters(ismember(parameterNames, 'Vfeces')).Value;
    Vsi = modelGutUniversal.Parameters(ismember(parameterNames, 'Vsi')).Value;
    VcecCVR = modelGutUniversal.Parameters(ismember(parameterNames, 'VcecumCVR')).Value;
    VcolCVR = modelGutUniversal.Parameters(ismember(parameterNames, 'VcolonCVR')).Value;
    VfecCVR = modelGutUniversal.Parameters(ismember(parameterNames, 'VfecesCVR')).Value;
    if ismember('Vurine', parameterNames)
       VurineDRUG = modelGutUniversal.Parameters(ismember(parameterNames, 'Vurine')).Value;
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % SIMULATE MODEL WITH data from datafile1 to estimate host parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %INITIAL PARAMETER INITIALIZATION 
    initPar = cell(nnz(parameterNames_fitdata==1),1);
    initParValue = zeros(nnz(parameterNames_fitdata==1),1);
    idx =1;
    for i=1:length(parameterNames)
        if parameterNames_fitdata(i)==1
           initPar{idx} = ['log(' parameterNames{i}, ')'];
           initParValue(idx) = parameterInitialValues(i);
           idx = idx+1;
        end
    end
    initParValue(initParValue==0) = rand(nnz(initParValue==0),1);
      
    % Possibly run the model several times to select the best fit
    nruns = 1; %number of runs to run the model
    randomRuns_init = zeros(nruns, length(initPar));
    randomRuns_results = cell(nruns,1);
    randomRunsBact_init = zeros(nruns, length(initPar));
    randomRunsBact_results = cell(nruns,1);

    for run_i = 1:nruns
        init = rand(1,length(initPar));
        for i=1:length(init)
            if initParValue(i)
                init(i) = initParValue(i);
            end
        end
        randomRuns_init(run_i,:) = init;
    end
    for run_i = 1:nruns
        init = randomRuns_init(run_i,:);
        estimated_parameters = estimatedInfo(initPar, 'InitialValue',init);%'
        responseMap = strcat(metNamesMap(useForFitting,1), ' = ',...
                         metNamesMap(useForFitting,2));

        % optimization step
        rng('default')
        globalMethod = 'ga';
        options = optimoptions(globalMethod);
        hybridMethod = 'fminsearch';
        hybridopts = optimset('Display', 'none');
        options = optimoptions(options, 'HybridFcn', {hybridMethod, hybridopts});
        resultsHost = sbiofit(modelGutUniversal, gd,responseMap,estimated_parameters,[],...
                          globalMethod,options,'pooled',true);
        % save fitting results
        randomRuns_results{run_i} = resultsHost;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set model parameters to optimized values and simulate
    resultsHost = randomRuns_results{run_i};
    for i=1:length(resultsHost.ParameterEstimates.Name)
        modelGutUniversal.Parameters(ismember(parameterNames, resultsHost.ParameterEstimates.Name(i))).Value =...
             resultsHost.ParameterEstimates.Estimate(i);
    end

    % define volumes according to the compartments in the model 
    if files_i<5
         curVolumes = [Vsi Vsi Vsi Vsi Vsi...
                  VbloodDRUG VbloodDRUG VcecGF VcolGF VfecGF...
                  VcecGF VcolGF VfecGF VbloodDRUG];
         
    else % adapt volumes for extended model
%         curVolumes = [Vsi Vsi Vsi Vsi Vsi...
%                   VbloodDRUG 1 VbloodDRUG VcecGF VcolGF VfecGF...
%                   Vsi Vsi Vsi...
%                   VcecGF VcolGF VfecGF 1 ... 
%                   VcecGF VcolGF VfecGF VbloodDRUG ...
%                   VcecGF VcolGF VfecGF VbloodDRUG...
%                   VurineDRUG ...
%                   VbloodDRUG 1 ...
%                   VurineDRUG VurineDRUG ...
%                   VurineDRUG 1 ...
%                   VbloodDRUG VurineDRUG ...
%                   VurineDRUG ...
%                   ];
    curVolumes = [Vsi Vsi ...
                  VbloodDRUG VurineDRUG ... 
                  Vsi Vsi Vsi...
                  VbloodDRUG VbloodDRUG...
                  1 1 ...
                  VcecGF VcolGF VfecGF ... 
                  VcecGF VcolGF VfecGF  ...
                  VbloodDRUG ...
                  VurineDRUG ...
                  VurineDRUG ...
                  VurineDRUG ...
                  Vsi Vsi Vsi ...
                  VbloodDRUG 1 1 ...
                  VcecGF VcolGF VfecGF ...
                  VcecGF VcolGF VfecGF...
                  VbloodDRUG...
                  VurineDRUG ...
                  VurineDRUG ...
                  ];
    end
 
              
    % simulate the model according to the set parameters
    simulate_pbpk_model(t, curVolumes, modelGutUniversal, ...
        metNamesMap, useForFitting)
    suptitle(datafilename1)
    
    % save model fitting and predictions to a file  
    prepare_data_for_tables_public(t, curVolumes, modelGutUniversal, resultsHost, ...
        metNamesMap, useForFitting, outfilename1, initParValue)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run the model with the data from datafile2 to estimate bacterial
    % coefficients
    modelGutUniversalBact = copyobj(modelGutUniversal);

    [t,metNamesMap, gd, useForFitting] = load_data_from_file(datafilename2);


    %INITIAL PARAMETER INITIALIZATION
    initPar = cell(nnz(parameterNames_fitdata==2),1);
    initParValue = 1+rand(nnz(parameterNames_fitdata==2),1);
    idx =1;
    for i=1:length(parameterNames)
        if parameterNames_fitdata(i)==2
           initPar{idx} = ['log(' parameterNames{i}, ')'];
           initParValue(idx) = parameterInitialValues(i);
           idx = idx+1;
        end
    end
    initParValue(initParValue==0) = rand(nnz(initParValue==0),1);

        estimated_parameters = estimatedInfo(initPar, 'InitialValue',initParValue);%'
        responseMap = strcat(metNamesMap(useForFitting,1), ' = ',...
                             metNamesMap(useForFitting,2));

        rng('default')
        globalMethod = 'ga';
        options = optimoptions(globalMethod);
        hybridMethod = 'fminsearch';
        hybridopts = optimset('Display', 'none');
        options = optimoptions(options, 'HybridFcn', {hybridMethod, hybridopts});
        resultsHostBact = sbiofit(modelGutUniversalBact, gd,responseMap,estimated_parameters,[],...
                          globalMethod,options,'pooled',true);

        randomRunsBact_init(run_i,1:length(init)) = init;
        randomRunsBact_results{run_i} = resultsHostBact;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set model parameters to optimized values and simulate
        parameterNames = cell(size(modelGutUniversalBact.Parameters));
        for i=1:length(modelGutUniversalBact.Parameters)
            parameterNames{i} = modelGutUniversalBact.Parameters(i).Name;
        end

        for i=1:length(resultsHostBact.ParameterEstimates.Name)
            modelGutUniversalBact.Parameters(ismember(parameterNames, resultsHostBact.ParameterEstimates.Name(i))).Value =...
                 resultsHostBact.ParameterEstimates.Estimate(i);
        end

        % define tissue weights depending on the model
        if files_i<5
            curVolumes = [Vsi Vsi Vsi Vsi Vsi...
                  VbloodDRUG VbloodDRUG VcecCVR VcolCVR VfecCVR...
                  VcecCVR VcolCVR VfecCVR VbloodDRUG];
        else
%             curVolumes = [Vsi Vsi Vsi Vsi Vsi...
%                   VbloodDRUG 1 VbloodDRUG VcecCVR VcolCVR VfecCVR...
%                   Vsi Vsi Vsi...
%                   VcecCVR VcolCVR VfecCVR 1 ... 
%                   VcecCVR VcolCVR VfecCVR VbloodDRUG ...
%                   VcecCVR VcolCVR VfecCVR VbloodDRUG...
%                   VurineDRUG ...
%                   VbloodDRUG 1 ...
%                   VurineDRUG VurineDRUG ...
%                   VurineDRUG 1 ...
%                   VbloodDRUG VurineDRUG ...
%                   VurineDRUG ...
%                   ];
                curVolumes = [Vsi Vsi ...
                  VbloodDRUG VurineDRUG ... 
                  Vsi Vsi Vsi...
                  VbloodDRUG VbloodDRUG...
                  1 1 ...
                  VcecCVR VcolCVR VfecCVR ... 
                  VcecCVR VcolCVR VfecCVR  ...
                  VbloodDRUG ...
                  VurineDRUG ...
                  VurineDRUG ...
                  VurineDRUG ...
                  Vsi Vsi Vsi ...
                  VbloodDRUG 1 1 ...
                  VcecCVR VcolCVR VfecCVR ...
                  VcecCVR VcolCVR VfecCVR...
                  VbloodDRUG...
                  VurineDRUG ...
                  VurineDRUG ...
                  ];
        end
        % simulate the model according to the set parameters
        simulate_pbpk_model(t, curVolumes, modelGutUniversalBact, ...
            metNamesMap, useForFitting)
        suptitle(datafilename2)
        % save model fitting and predictions to a file  
        prepare_data_for_tables_public(t, curVolumes, modelGutUniversalBact, resultsHostBact, ...
            metNamesMap, useForFitting, outfilename2, initParValue)
        if perform_sensitivity_analysis_flag
            perform_global_sensitivity_analysis(modelGutUniversal,...
                                             resultsHost,...
                                             resultsHostBact,...
                                             1000,...
                                             0,... % change F flag
                                             0,... % change dose flag
                                             0); % simulate within estimated flag
        end
end