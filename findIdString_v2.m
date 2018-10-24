%USAGE
%findIdString( model.rxns, 'biomass_synthesis')
%findIdString( model.met, 'AKGc')
%findIdString(model, 'string') search 'string' in all fields


function [ Id_string ] = findIdString_v2( where, string)
% funzione che individua l'ID della reazione o del metabolita a partire dal nome

if ~isstruct(where)
    Id_string = find(strcmp(where, string)==1);
    if isempty(Id_string)
        Id_string = 0;
    end
else
    Fields = [];
    SNames = fieldnames(where);
    for loopIndex = 1:numel(SNames)
        Id_string = find(strcmp(where.(SNames{loopIndex}), string)==1);
        if ~isempty(Id_string)
            Fields = [Fields SNames{loopIndex} num2cell(Id_string')];
        end
    end
    if ~isempty(Fields)
    Id_string = Fields;
    else
        Id_string = 0;
    end
end
end

