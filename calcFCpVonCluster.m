function [tableOut] = calcFCpVonCluster(structCluster, varargin)
% Compute pValue with many parametric or not parametric tests and the Fold Change 
% of the mean of the transcript between the samples in the two or more
% cluster.
% If there is more than two cluster the comparison on each reactions could be performed
% between all the combination from the n cluster (method 'combinatorial') or between each cluster
% and the mean of the (mean) values from all cluster (method 'single').
%
%
% USAGE:
%
%    [tableOut] = calcFCpVonCluster(structCluster, 'method', 'combinatorial', 'test', 'ttest2', 'param', [{'tail'},{'both'}], 'normalization', 'max', 'fc', true, 'log2fc', true, 'absFC', true, 'header', true, 'filename', 'TableOut.txt');
%
% INPUT:
%    structCluster      Struct with many fields as clusters created with
%                       extractCluster function.
%
% OPTIONAL INPUTS:
%    method:            'single', 'combinatorial' (Default = 'combinatorial')
%    test:              Statistical test to perform:  
%                       'anova1','kruskalwallis','kstest2','ttest2','vartest2', 'vartestn', 'zscore', 'none'
%                       (default = ttest2)
%    param:             Parameters pass to some functions that need it i.e. ttest2 (Default [])
%    normalization:     Normalization of the mean value of each reactions in each clusters 
%                       'max': divide each values with the max of the mean
%                       values
%                       'sum': divide each values with the sum of the value means
%                       (default = 'none'): do nothing
%   fc:                 compute the fold change if true, nothing if is false
%   log2fc:             if true and if fc = true log2 transform the fc output
%   header:             true if the first row contain the samples id.
%                       Please note that the first column is always use as
%                       reactions id.
%   filename:           if passed the output table will be saved with that
%                       filename in a table separated file (Default: '')
%
% OUTPUTS:
%    tableOut:         	Table with the resoults.
%

% parser input
p = inputParser;
p.CaseSensitive = false;
defaultMethod = 'combinatorial';
expectedMehod = {'single', 'combinatorial'};
defaultTest = 'ttest2';
expectedTest = {'anova1','kruskalwallis','kstest2', 'ttest2','vartest2', 'vartestn', 'zscore', 'none'};
defaultParam = []; 
defaultFc =  true;
defaultLog2FC = true;
defaultAbsFC = false;
defaultNormalization = 'none';
expectedNormalization = {'max','sum','none'};
defaultHeader =  true;
defaultFileName = '';

addRequired(p,'structCluster',@isstruct);
addParameter(p,'method',defaultMethod,@(x) any(validatestring(x,expectedMehod)));
addParameter(p,'test',defaultTest,@(x) any(validatestring(x,expectedTest)));
addParameter(p,'param',defaultParam, @iscell);
addParameter(p,'fc', defaultFc, @islogical);
addParameter(p,'log2FC',defaultLog2FC,@islogical);
addParameter(p,'absFC',defaultAbsFC,@islogical);
addParameter(p,'normalization',defaultNormalization,@(x) any(validatestring(x,expectedNormalization)));
addParameter(p,'header',defaultHeader,@islogical);
addParameter(p, 'filename', defaultFileName, @ischar);

parse(p,structCluster,varargin{:});

structCluster = p.Results.structCluster;
method = p.Results.method;
test = p.Results.test;
param = p.Results.param;
FCflag =  p.Results.fc;
log2FC = p.Results.log2FC;
absFC = p.Results.absFC;
normalization = p.Results.normalization;
header = p.Results.header;
filename = p.Results.filename;

tableOut = table();
ClusterName = fieldnames(structCluster);
if header
    tableOut.('VarName') = structCluster.(ClusterName{1})(2:end,1);
else
    tableOut.('VarName') = structCluster.(ClusterName{1})(1:end,1);
end
for i=1:length(ClusterName)
    %first column with reaction names remove
    structCluster.(ClusterName{i})(:,1)=[];
    if header %first row with sample names remove
        structCluster.(ClusterName{i})(1,:)=[];
    end
    structCluster.(ClusterName{i}) = cell2mat(structCluster.(ClusterName{i}));
end

if strcmp(method, 'single')
    tableOut = [tableOut, methodSingle(structCluster, test, param, FCflag, log2FC, absFC, normalization)];
elseif (strcmp(method, 'combinatorial'))
    tableOut = [tableOut, methodComb(structCluster, test, param, FCflag, log2FC, absFC, normalization)];
else
    warning('Invalid method: ''single'' or ''combinatorial'' allowed');
    return
end

if ~isempty(filename)
    writetable(tableOut, filename, 'delimiter', '\t');
end

end


function [tableOut] = methodSingle(structCluster, test, param, FCflag, log2FC, absFC, normalization)
tableOut = table();
ClusterName = fieldnames(structCluster);
nClust =  length(ClusterName);
nRxn = size(structCluster.(ClusterName{1}), 1);

if ~strcmp(test, 'none')
    switch  test
        case 'anova1'
            %performs a one-way ANOVA for comparing the means of two or more
            %groups of data. It returns the p-value for the null hypothesis that the
            %means of the groups are equal.
            pV = Anova1onCluster(structCluster);
            
        case 'kruskalwallis'
            %Kruskal-Wallis test. Tests if multiple samples are all drawn from the
            %same populations (or equivalently, from different populations with the
            %same distribution), against the alternative that they are not all drawn
            %from the same population.
            pV = KWonCluster(structCluster);
            
        case 'kstest2'
            %Two-sample Kolmogorov-Smirnov test. Tests if two samples come from the
            %same continuous distribution, against the alternative that they do not
            %come from the same distribution.
            if nClust ~= 2
                error('Kstest2 compare two samples, please use ''combinatorial'' method or kruskalwallis instead');
            else
                pV = KSonCluster(structCluster, param);
            end
            
            
        case 'ttest2'
            %Two-sample t-test. Tests if two independent samples come from normal
            %distributions with unknown but equal (or, optionally, unequal) variances
            %and the same mean, against the alternative that the means are unequal.
            if nClust ~= 2
                error('ttest2 compare two samples, please use ''combinatorial'' method or a multiple samples test instead');
            else
                pV = TSonCluster(structCluster, param);
            end
            
        case 'vartest2'
            %Two-sample F-test for equal variances. Tests if two independent samples
            %come from normal distributions with the same variance,
            %against the alternative that they come from normal distributions with
            %different variances.
            if nClust ~= 2
                error('vartest2 compare two samples, please use ''combinatorial'' method or vartestn instead');
            else
                pV = VT2onCluster(structCluster, param);
            end
            
        case 'vartestn'
            %Bartlett multiple-sample test for equal variances. Tests if multiple
            %samples come from normal distributions with the same variance, against
            %the alternative that they come from normal distributions with different
            %variances.
            pV = VTNonCluster(structCluster, param);
            
        case 'zscore'
            %Implementato anche lo Z Score come nella funzione getZscores_v2
            % i param sono woErr o wErr per calcolarlo senza o con
            % Errore Standard
            error('Z Score availabel only with ''combinatorial'' method');
            
    end
    tableOut.(test) = pV;
end

% fold change compute respect to the mean of the values

if FCflag || ~strcmp(normalization, 'none')
    meanC_i = zeros(nRxn, nClust);
    for i=1:nClust
        meanC_i(:,i) = nanmean(structCluster.(ClusterName{i}),2);
    end
    meanC_T = nanmean(meanC_i, 2);
    maxC_T = nanmax(meanC_i, [], 2);
    
    if ~strcmp(normalization, 'none')
        switch normalization
            case 'max' 
                norm_i = meanC_i ./ maxC_T;
            case 'sum' 
                norm_i = meanC_i ./nansum(meanC_i,2);
        end
        for i=1:nClust
            tableOut.(strcat('Norm_', ClusterName{i})) = norm_i(:,i);
        end
    end
    
    if FCflag
        if log2FC
            FC_i = log2(meanC_i) - log2(meanC_T);
        else
            FC_i = meanC_i ./ meanC_T;
        end
        % save the resoult in a table
        
        for i=1:nClust
            tableOut.(strcat('FC_', ClusterName{i})) = FC_i(:,i);
            if absFC
                tableOut.(strcat('absFC_', ClusterName{i})) = abs(FC_i(:,i));
            end
        end
    end
end

end

function [tableOut] = methodComb(structCluster, test, param, FCflag, log2FC, absFC, normalization)
tableOut = table();
ClusterName = fieldnames(structCluster);
ClustComb = combnk(ClusterName,2);
nComb = size(ClustComb, 1);
nRxn = size(structCluster.(ClusterName{1}), 1);
for k=1:nComb
    structTmp.('A') = structCluster.(ClustComb{k,1});
    structTmp.('B') = structCluster.(ClustComb{k,2});
    nameA = ClustComb{k,1};
    nameB = ClustComb{k,2};
    LabelTableOut = strcat(nameA, '_vs_', nameB);
    nClust = 2;
    
    if ~strcmp(test, 'none')
        switch  test
            case 'anova1'
                %performs a one-way ANOVA for comparing the means of two or more
                %groups of data. It returns the p-value for the null hypothesis that the
                %means of the groups are equal.
                pV = Anova1onCluster(structTmp);
                
            case 'kruskalwallis'
                %Kruskal-Wallis test. Tests if multiple samples are all drawn from the
                %same populations (or equivalently, from different populations with the
                %same distribution), against the alternative that they are not all drawn
                %from the same population.
                warning('With ''combinatorial'' method we suggest to use kstest2 instead');
                pV = KWonCluster(structTmp);
                
            case 'kstest2'
                %Two-sample Kolmogorov-Smirnov test. Tests if two samples come from the
                %same continuous distribution, against the alternative that they do not
                %come from the same distribution.
                if nClust ~= 2
                    error('Kstest2 compare two samples, please use ''combinatorial'' method or kruskalwallis instead');
                else
                    pV = KSonCluster(structTmp, param);
                end
                
            case 'ttest2'
                %Two-sample t-test. Tests if two independent samples come from normal
                %distributions with unknown but equal (or, optionally, unequal) variances
                %and the same mean, against the alternative that the means are unequal.
                if nClust ~= 2
                    error('ttest2 compare two samples, please use ''combinatorial'' method or a multiple samples test instead');
                else
                    pV = TSonCluster(structTmp, param);
                end
                
            case 'vartest2'
                %Two-sample F-test for equal variances. Tests if two independent samples
                %come from normal distributions with the same variance,
                %against the alternative that they come from normal distributions with
                %different variances.
                if nClust ~= 2
                    error('vartest2 compare two samples, please use ''combinatorial'' method or vartestn instead');
                else
                    pV = VT2onCluster(structTmp, param);
                end
                
            case 'vartestn'
                %Bartlett multiple-sample test for equal variances. Tests if multiple
                %samples come from normal distributions with the same variance, against
                %the alternative that they come from normal distributions with different
                %variances.
                warning('With ''combinatorial'' method we suggest to use vartest2 instead');
                pV = VTNonCluster(structTmp, param);
                
            case 'zscore'
                %Implementato anche lo Z Score come nella funzione getZscores_v2
                % i param sono woErr o wErr per calcolarlo senza o con
                % Errore Standard
                pV = pvZscore(structTmp, param);
        end
        tableOut.(strcat(test, '_pV_', LabelTableOut)) = pV;
    end
    
    if FCflag || ~strcmp(normalization, 'none')
        meanC = zeros(nRxn, 2);
        % compute fold change between A and B
        meanC(:,1) = nanmean(structCluster.(nameA),2);
        meanC(:,2) = nanmean(structCluster.(nameB),2);
    end
    if ~strcmp(normalization, 'none')
        switch normalization
            case 'max' 
                norm_T = meanC ./ nanmax(meanC, [], 2);
            case 'sum' 
                norm_T = meanC ./nansum(meanC,2);
        end
        
        tableOut.(strcat('Norm_',normalization, '_', LabelTableOut)) = norm_T;
    end
    
    if FCflag
        if log2FC
            meanC = log2(meanC);
            FC = meanC(:,1) - meanC(:,2);
        else
            FC = meanC(:,1)./meanC(:,2);
        end
        tableOut.(strcat('FC_',LabelTableOut)) = FC;
        if absFC
            tableOut.(strcat('absFC_',LabelTableOut)) = abs(FC);
        end
    end
end
end


function [pV] = KWonCluster(structCluster)
ClusterName = fieldnames(structCluster);
nClust = length(ClusterName);
nRxn = size(structCluster.(ClusterName{1}), 1);
pV = ones(nRxn, 1);
for j=1:nRxn
    Aggregate = [];
    ClustLab = [];
    for i=1:nClust
        tmpAggregate = structCluster.(ClusterName{i})(j,:)';
        tmpClustLab = zeros(size(structCluster.(ClusterName{i}),2),1);
        tmpClustLab(:) = i;
        Aggregate = [Aggregate; tmpAggregate];
        ClustLab = [ClustLab; tmpClustLab];
    end
    pV(j) = kruskalwallis(Aggregate, ClustLab, 'off');
end
end

function [pV] = Anova1onCluster(structCluster)
ClusterName = fieldnames(structCluster);
nClust = length(ClusterName);
nRxn = size(structCluster.(ClusterName{1}), 1);
pV = ones(nRxn, 1);
for j=1:nRxn
    Aggregate = [];
    ClustLab = [];
    for i=1:nClust
        tmpAggregate = structCluster.(ClusterName{i})(j,:)';
        tmpClustLab = zeros(size(structCluster.(ClusterName{i}),2),1);
        tmpClustLab(:) = i;
        Aggregate = [Aggregate; tmpAggregate];
        ClustLab = [ClustLab; tmpClustLab];
    end
    pV(j) = anova1(Aggregate, ClustLab, 'off');
end
end

function [pV] = KSonCluster(structCluster, param)
ClusterName = fieldnames(structCluster);
nRxn = size(structCluster.(ClusterName{1}), 1);
pV = ones(nRxn, 1);
for j=1:nRxn
    try
        if isempty(param)
            [~,pV(j)] = kstest2(structCluster.(ClusterName{1})(j,:), structCluster.(ClusterName{2})(j,:));
        else
            [~,pV(j)] = kstest2(structCluster.(ClusterName{1})(j,:), structCluster.(ClusterName{2})(j,:), param{:});
        end
    catch ME
        if(strcmp(ME.identifier, 'stats:kstest2:NotEnoughData'))
            pV(j) = NaN;
        else
            error('kstet2, wrong param input');
        end
    end
end
end

function [pV] = TSonCluster(structCluster, param)
ClusterName = fieldnames(structCluster);
nRxn = size(structCluster.(ClusterName{1}), 1);
pV = ones(nRxn, 1);
for j=1:nRxn
    % verify that all the value are not the same. If so the pValue is 1.
    % some test return a error if all the value are the same
    if range([structCluster.(ClusterName{2})(j,:), structCluster.(ClusterName{1})(j,:)]) == 0
        pV(j) = 1;
    else
        try
            if isempty(param)
                [~,pV(j)] = ttest2(structCluster.(ClusterName{1})(j,:), structCluster.(ClusterName{2})(j,:));
            else
                [~,pV(j)] = ttest2(structCluster.(ClusterName{1})(j,:), structCluster.(ClusterName{2})(j,:), param{:});
            end
        catch ME
            if(strcmp(ME.identifier, 'stats:ttest2:NotEnoughData'))
                pV(j) = NaN;
            else
                error('ttest2, wrong param input');
            end
        end
    end
end
end

function [pV] = VT2onCluster(structCluster, param)
ClusterName = fieldnames(structCluster);
nRxn = size(structCluster.(ClusterName{1}), 1);
pV = ones(nRxn, 1);
for j=1:nRxn
    % verify that all the value are not the same. If so the pValue is 1.
    % some test return a error if all the value are the same
    if range([structCluster.(ClusterName{2})(j,:), structCluster.(ClusterName{1})(j,:)]) == 0
        pV(j) = 1;
    else
        try
            if isempty(param)
                [~,pV(j)] = vartest2(structCluster.(ClusterName{1})(j,:), structCluster.(ClusterName{2})(j,:));
            else
                [~,pV(j)] = vartest2(structCluster.(ClusterName{1})(j,:), structCluster.(ClusterName{2})(j,:), param{:});
            end
        catch ME
            if(strcmp(ME.identifier, 'stats:vartest2:NotEnoughData'))
                pV(j) = NaN;
            else
                error('vartest2, wrong param input');
            end
        end
    end
end
end

function [pV] = VTNonCluster(structCluster, param)
param = [{'display'},{'off'}, param];
ClusterName = fieldnames(structCluster);
nClust = length(ClusterName);
nRxn = size(structCluster.(ClusterName{1}), 1);
pV = ones(nRxn, 1);
for j=1:nRxn
    Aggregate = [];
    ClustLab = [];
    for i=1:nClust
        tmpAggregate = structCluster.(ClusterName{i})(j,:)';
        tmpClustLab = zeros(size(structCluster.(ClusterName{i}),2),1);
        tmpClustLab(:) = i;
        Aggregate = [Aggregate; tmpAggregate];
        ClustLab = [ClustLab; tmpClustLab];
    end
    % verify that all the value are not the same. If so the pValue is 1.
    % some test return a error if all the value are the same
    if range(Aggregate) == 0
        pV(j) = 1;
    else
        pV(j) = vartestn(Aggregate, ClustLab, param{:});
    end
end
end

function [pV] = pvZscore(structCluster, param)
ClusterName = fieldnames(structCluster);
nRxn = size(structCluster.(ClusterName{1}), 1);
pV = ones(nRxn, 1);
for j=1:nRxn
    v1 = structCluster.(ClusterName{1})(j,:);
    v2 = structCluster.(ClusterName{2})(j,:);
    
    % verify that all the value are not the same. If so the pValue is 1.
    % some test return a error if all the value are the same
    if range([v1,v2]) == 0
        pV(j) = 1;
    else
        
        if isempty(param) || strcmp(param, 'woErr')
            Zi=( mean(v2) - mean(v1) ) ./ sqrt(var(v2)+var(v1));
        elseif strcmp(param, 'wErr')
            Zi=( mean(v2) - mean(v1) ) ./ sqrt(var(v2)./length(v2) + var(v1)./length(v1));
        else
            error('Zscore, wrong param input');
        end
        
        pV(j) = 1 - normcdf(abs(Zi));
    end
end
end





