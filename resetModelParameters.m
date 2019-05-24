function newModel = resetModelParameters(baseModel, Fdrug, dose)
% reset model parameters to 0
% set species initial values according to Fdrug and dose
    newModel = copyobj(baseModel);
    % reset initial P_si and P_si0 concentrations 
    % to the provided Fdrug and dise values
    P_si = sbioselect(newModel, 'Name','P_si');
    P_si.InitialAmount = Fdrug*dose;
    P_si0 = sbioselect(newModel, 'Name','P_si0');
    P_si0.InitialAmount = (1-Fdrug)*dose;
    % reset values of all changing parameters to 0 (Notes~=0)
    for i=1:length(newModel.Parameters)
        if ~isequal(newModel.Parameters(i).Tag, 'constant')
            newModel.Parameters(i).Value = 0;
        end
    end