%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simulate_pbpk_model(t, curVolumes, modelGutUniversal, ...
    metNamesMap, useForFitting)

%% Simulate
    configsetObj = getconfigset(modelGutUniversal);
    configsetObj.SolverOptions.SensitivityAnalysis = false;
    set(configsetObj, 'StopTime', max(t.Time)+1);
    simData = sbiosimulate(modelGutUniversal);
   
    fig = figure('units','normalized','outerposition',[0 0 1 1]);
    lineWidth = 0.5;
    for i=1:length(simData.dataNames)
        subplot(3,ceil(length(simData.dataNames)/3),i)
        hold on
        if ismember(simData.dataNames{i}, metNamesMap(useForFitting,1)) % used for fitting
            h1 = plot(simData.Time,simData.Data(:,i)/curVolumes(i), 'LineWidth', 2);
        else
            h2 = plot(simData.Time,simData.Data(:,i)/curVolumes(i), '--', 'LineWidth', 2);
        end
        title(simData.DataNames{i})
        xlim([0 10])
        axis square
    end
    mycolors = [0 .6 .5];  % bluish green    --- Bth WT
    %mycolors = [0 .45 .7]; % blue            --- mouse /CVR
    %plot experimental data on top
    for i=1:length(simData.dataNames)
        subplot(3,ceil(length(simData.dataNames)/3),i)
        hold on
        if nnz(ismember(metNamesMap(:,1), simData.dataNames{i}))
            dataMet = metNamesMap(ismember(metNamesMap(:,1), simData.dataNames{i}),2);
               for j=1:length(t.Time)
                    if ~isnan(t{j, dataMet})
                        h3 = plot([t.Time(j)-lineWidth, t.Time(j)+lineWidth],...
                                     [t{j, dataMet}, t{j, dataMet}]/curVolumes(i), 'LineWidth', 3, 'Color', mycolors);
                    end
                end
        end
    end
    % plot serum metabolite if bacterial contribution is nonzero
    if nnz(simData.Data(:,ismember(simData.dataNames, 'M_serumbact'))>0)
        subplot(3,ceil(length(simData.dataNames)/3),find(ismember(simData.dataNames, 'M_serum')))
        hold on
        plot(simData.Time, simData.Data(:,ismember(simData.dataNames, 'M_serum'))/...
                                        curVolumes(ismember(simData.dataNames, 'M_serum')) + ...
                           simData.Data(:,ismember(simData.dataNames, 'M_serumbact'))/...
                                        curVolumes(ismember(simData.dataNames, 'M_serum')), 'r', 'LineWidth', 2)
        plot(simData.Time, simData.Data(:,ismember(simData.dataNames, 'M_serumbact'))/...
                                        curVolumes(ismember(simData.dataNames, 'M_serum')), '--r', 'LineWidth', 2)
    end
    if nnz(simData.Data(:,ismember(simData.dataNames, 'M1_serumbact'))>0)
        subplot(3,ceil(length(simData.dataNames)/3),find(ismember(simData.dataNames, 'M1_serum')))
        hold on
        plot(simData.Time, simData.Data(:,ismember(simData.dataNames, 'M1_serum'))/...
                                        curVolumes(ismember(simData.dataNames, 'M1_serum')) + ...
                           simData.Data(:,ismember(simData.dataNames, 'M1_serumbact'))/...
                                        curVolumes(ismember(simData.dataNames, 'M1_serum')), 'r', 'LineWidth', 2)
        plot(simData.Time, simData.Data(:,ismember(simData.dataNames, 'M1_serumbact'))/...
                                        curVolumes(ismember(simData.dataNames, 'M1_serum')), '--r', 'LineWidth', 2)
    end
    if nnz(simData.Data(:,ismember(simData.dataNames, 'M2_serumbact'))>0)
        subplot(3,ceil(length(simData.dataNames)/3),find(ismember(simData.dataNames, 'M2_serum')))
        hold on
        plot(simData.Time, simData.Data(:,ismember(simData.dataNames, 'M2_serum'))/...
                                        curVolumes(ismember(simData.dataNames, 'M2_serum')) + ...
                           simData.Data(:,ismember(simData.dataNames, 'M2_serumbact'))/...
                                        curVolumes(ismember(simData.dataNames, 'M2_serum')), 'r', 'LineWidth', 2)
        plot(simData.Time, simData.Data(:,ismember(simData.dataNames, 'M2_serumbact'))/...
                                        curVolumes(ismember(simData.dataNames, 'M2_serum')), '--r', 'LineWidth', 2)
    end
    if nnz(simData.Data(:,ismember(simData.dataNames, 'M3_serumbact'))>0)
        subplot(3,ceil(length(simData.dataNames)/3),find(ismember(simData.dataNames, 'M3_serum')))
        hold on
        plot(simData.Time, simData.Data(:,ismember(simData.dataNames, 'M3_serum'))/...
                                        curVolumes(ismember(simData.dataNames, 'M3_serum')) + ...
                           simData.Data(:,ismember(simData.dataNames, 'M3_serumbact'))/...
                                        curVolumes(ismember(simData.dataNames, 'M3_serum')), 'r', 'LineWidth', 2)
        plot(simData.Time, simData.Data(:,ismember(simData.dataNames, 'M3_serumbact'))/...
                                        curVolumes(ismember(simData.dataNames, 'M3_serum')), '--r', 'LineWidth', 2)
    end
    set(groot, 'DefaultTextInterpreter', 'none')
    legend([h1 h2 h3], {'Model fit', 'Prediction', 'Experimental data'})