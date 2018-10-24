function [ClustersOut] = extractCluster(dataIn, matchList)
%
% USAGE:
%
%   [ClustersOut] = extractCluster(dataIn, matchList)
%
% INPUT:
%   dataIn:         matlab table or cell array with each cols contain an osservation except
%                   for the first Col which contain variable names. If is a
%                   cell array first row contain samples id to match
%                   matchList.
%   matchList:      Cell array with osservation to extract. If there is two
%                   column the second will use as cluster labels.
%                   Please note that this labels can not be number but only
%                   string.
%
% OUTPUTS:
%   ClustersOut:    Dataset with only columns in match list or a struct
%                   with many fields as different labels passed.

ClustersOut = struct();
nCols = size(matchList, 2);
if istable(dataIn)
    rowsName = cellstr(dataIn.(1));
    datasetIn = dataIn{:,2:end}';
    varNames = dataIn.Properties.VariableNames;
else
    rowsName = cellstr(dataIn(2:end,1));
    datasetIn = cell2mat(dataIn(2:end,2:end))';
    varNames = cellstr(dataIn(1,:));
end
if nCols > 1
    %extract label in col 2.
    dataLabel = unique(matchList(:,2));
    nClass = length(dataLabel);
    clustL = 1;
    for i=1:nClass
        idxL = ismember(matchList(:,2), dataLabel(i));
        ListaItem = matchList(idxL,1);
        idxI = ismember(varNames, ListaItem);
        idxI(1) = 1; %aggiunge la corrispondenza sulla prima colonna (quella che contiene le reazioni)
        extracted = datasetIn(idxI(2:end)==1,:)';
        extracted = [rowsName num2cell(extracted)];
        extracted = [varNames(idxI); extracted];
        %dataLabel{i}(~ismember(double(dataLabel{i}),[65:90 97:122])) = ''; %rimuove caratteri speciali dalla cluster label
        try
        ClustersOut.(dataLabel{i}) = extracted;
        catch ME
            if (strcmp(ME.identifier,'MATLAB:AddField:InvalidFieldName'))
                ClustersOut.(num2str(clustL, 'clust%i')) =  extracted;
                disp(['Invalid label: ' dataLabel{i} '. Replaced with: ' num2str(clustL, 'clust%i')])
                clustL = clustL + 1;
            else
                rethrow(ME);
            end
        end
    end
else
    idxI = ismember(varNames, matchList);
    idxI(1) = 1; %aggiunge la corrispondenza sulla prima colonna (quella che contiene le reazioni)
    extracted = datasetIn(idxI(2:end)==1,:)';
    extracted = [rowsName num2cell(extracted)];
    extracted = [varNames(idxI==1); extracted];
    ClustersOut = extracted;
end