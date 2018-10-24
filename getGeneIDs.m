% la funzione contatta genenames.org e chiede per ogni gene tutte le
% informazioni sui vari ID disponibili.

function [tableComplete] = getGeneIDs(geneList, fromFeat, toFeat, fileName)

% geneList =  cell array [{'xxx'}, ] con i geni da cercare, il formato può
% essere una feature tra queste: {'vega_id'},{'locus_group'}, {'alias_symbol'}, {'prev_name'}, {'refseq_accession'}, {'hgnc_id'}, {'entrez_id'}, {'symbol'}, {'name'}, {'mgd_id'}, {'prev_symbol'}, {'alias_name'}, {'status'}, {'locus_type'}, {'rgd_id'}, {'rna_central_ids'}, {'ensembl_gene_id'}, {'ucsc_id'}, {'uniprot_ids'}, {'ena'}, {'ccds_id'}
% fromFeat =  cell array [{'xxx'}, ] contenente le feature con cui cercare nel database.
% toFeat = cell array [{'xxx'}, ] contenente le feature da salvare.
% fileName = se si vuole salvare l'output in un file tsv inserire il nome 'xxx.txt'


availableFeat = [{'hgnc_id'},{'symbol'},{'name'},{'status'},{'locus_type'},{'prev_symbol'},{'prev_name'},{'alias_symbol'},{'alias_name'},{'location'},{'date_approved_reserved'},{'date_modified'},{'date_name_changed'},{'ena'},{'entrez_id'},{'mgd_id'},{'iuphar'},{'cosmic'},{'orphanet'},{'pubmed_id'},{'refseq_accession'},{'gene_family'},{'date_symbol_changed'},{'vega_id'},{'lsdb'},{'ensembl_gene_id'},{'ccds_id'},{'locus_group'},{'omim_id'},{'uniprot_ids'},{'ucsc_id'},{'rgd_id'},{'gene_family_id'},{'location_sortable'},{'uuid'},{'x_version_'}];


if nargin < 3
    toFeat = [{'symbol'}, {'hgnc_id'}, {'ensembl_gene_id'}, {'entrez_id'}, {'uniprot_ids'}, {'name'}];
end
if strcmp(fromFeat, 'all') % cerca in tutti i campi
    fromFeat = availableFeat;
end


idxWrFeat = find(ismember(toFeat, availableFeat)==0);

if ~isempty(idxWrFeat)
    warning('Following feature(s) not available. Will be ignored');
    disp(toFeat(idxWrFeat));
    toFeat(idxWrFeat) = [];
    if isempty(toFeat)
        error('All features are not available');
    end
end
idxWrFeat = find(ismember(fromFeat, availableFeat)==0);
fromFeat(idxWrFeat) = [];

stQuery = 'http://rest.genenames.org/fetch/';
GenToProcess = 1:length(geneList);
genSt = 1;
f = 1;
errorQ = false;
tableComplete = table();
GenDone = [];
tmpFileName = cell('');

if isnumeric(geneList)
    entrez = true;
    geneList = num2cell(geneList);
else
    entrez = false;
end
h = waitbar(0,'Gene converted: 0%');
while ~isempty(GenToProcess)
    tableRes = table();
    for i=1:length(GenToProcess)
        if nargin < 2 || isempty(fromFeat)
            if entrez
                fromFeat = {'entrez_id'};
            elseif strncmp(geneList{GenToProcess(i)}, 'HGNC', 4)
                fromFeat = {'hgnc_id'};
            elseif strncmp(geneList{GenToProcess(i)}, 'ENSG', 4)
                fromFeat = {'ensembl_gene_id'};
            else
                fromFeat = [{'symbol'}, {'prev_symbol'},{'alias_symbol'}];
            end
        end
        fm = length(fromFeat);
        if ismember(fromFeat{f}, [{'entrez_id'}, {'kznf_gene_catalog'}, {'omim_id'}])
            query = strcat(stQuery, fromFeat{f}, '/', num2str(geneList{GenToProcess(i)}));
        else
            query = strcat(stQuery, fromFeat{f}, '/', geneList{GenToProcess(i)});
        end
        try
            data = webread(query);
        catch
            %check internet connection
            tf =  false;
            while tf == false
                try
                    address = java.net.InetAddress.getByName('www.google.com');
                    tf = true;
                catch
                    disp('No internet connection');
                    java.lang.Thread.sleep(1000) % if fail wait for one second
                end
            end
            if ~isempty(tableRes)
                tmpFileName = [tmpFileName; strcat('RisParz_From_gen_',num2str(genSt), '_to_', num2str(GenToProcess(i-1)),'.txt')];
                writetable(tableRes,tmpFileName{end},'Delimiter','\t');
                tableComplete = [tableComplete; tableRes];
            end
            if length(GenToProcess) == 1
                genSt = 1;
            else
                if length(GenToProcess) == 1
                    genSt = 1;
                else
                    genSt = GenToProcess(i+1);
                end
            end
            tableRes = table();
            continue
        end
        while data.response.numFound == 0 && f < fm
            f = f + 1;
            if ismember(fromFeat{f}, [{'entrez_id'}, {'kznf_gene_catalog'}, {'omim_id'}])
                query = strcat(stQuery, fromFeat{f}, '/', num2str(geneList{GenToProcess(i)}));
            else
                query = strcat(stQuery, fromFeat{f}, '/', geneList{GenToProcess(i)});
            end
            try
                data = webread(query);
            catch
                %check internet connection
                tf =  false;
                while tf == false
                    try
                        address = java.net.InetAddress.getByName('www.google.com');
                        tf = true;
                    catch
                        disp('No internet connection');
                        java.lang.Thread.sleep(1000) % if fail wait for one second
                    end
                end
                if ~isempty(tableRes)
                    tmpFileName = [tmpFileName; strcat('RisParz_From_gen_', num2str(genSt), '_to_', num2str(GenToProcess(i-1)),'.txt')];
                    writetable(tableRes,tmpFileName{end},'Delimiter','\t');
                    tableComplete = [tableComplete; tableRes];
                end
                errorQ =  true;
                break
            end
        end
        f = 1;
        if errorQ
            if length(GenToProcess) == 1
                genSt = 1;
            else
                genSt = GenToProcess(i+1);
            end
            tableRes = table();
            errorQ = false;
            continue
        end
        GenDone = [GenDone; GenToProcess(i)];
        %         if data.response.numFound > 1
        %             PrintFound = ['ambiguous: ' num2str(data.response.numFound) ' found'];
        %             symbol = {PrintFound};
        %             hgnc_id = {'ambiguous'};
        %             ensembl_gene_id = {'ambiguous'};
        %             entrez_id = {'ambiguos'};
        %             uniprot_ids = {'ambiguous'};
        %         else
        %             if isfield(data.response.docs, 'symbol')
        %                 symbol = {data.response.docs.symbol};
        %             else
        %                 symbol = {'not found'};
        %             end
        %             if isfield(data.response.docs, 'hgnc_id')
        %                 hgnc_id = {data.response.docs.hgnc_id};
        %             else
        %                 hgnc_id = {'not found'};
        %             end
        %             if isfield(data.response.docs, 'ensembl_gene_id')
        %                 ensembl_gene_id = {data.response.docs.ensembl_gene_id};
        %             else
        %                 ensembl_gene_id = {'not found'};
        %             end
        %             if isfield(data.response.docs, 'entrez_id')
        %                 entrez_id = {data.response.docs.entrez_id};
        %             else
        %                 entrez_id = {'not found'};
        %             end
        %             if isfield(data.response.docs, 'uniprot_ids')
        %                 uniprot_ids = data.response.docs.uniprot_ids;
        %             else
        %                 uniprot_ids = {'not found'};
        %             end
        %         end
        %         tabletmp = table(GenToProcess(i), geneList(GenToProcess(i)), symbol(1), hgnc_id(1), ensembl_gene_id(1), entrez_id(1), uniprot_ids(1),  'VariableNames', {'idx', 'Input', 'Symbol', 'HGNC', 'ENS', 'Entrez', 'Uniprot'});
        numGenPro = GenToProcess(i)/length(geneList);
        waitbar(numGenPro, h, ['Gene converted: ' num2str(round(numGenPro*100, 2)) '%']);
        tabletmp = table(GenToProcess(i), geneList(GenToProcess(i)),  'VariableNames', {'idx', 'Input'});
        tabletmp = [tabletmp, extractFeature(data.response, toFeat)];
        tableRes = [tableRes; tabletmp];
    end
    tmpFileName = [tmpFileName; strcat('RisParz_From_gen_',num2str(genSt), '_to_', num2str(GenToProcess(i)),'.txt')];
    writetable(tableRes,tmpFileName{end},'Delimiter','\t')
    GenToProcess = setdiff(GenToProcess, GenDone);
    tableComplete = [tableComplete; tableRes];
end
tableComplete = sortrows(tableComplete, 1);
if nargin >= 4
    writetable(tableComplete, fileName, 'Delimiter','\t');
end
%eliminare i file temporanei
delete(tmpFileName{:});

close(h);
end

function [tableTmp] = extractFeature(response, toFeat)

tableTmp = table();

if response.numFound > 1
    for i=1:length(toFeat)
        tableTmp.(toFeat{i}) = strcat('ambiguous: ', num2str(response.numFound), ' found');
    end
else
    for i=1:length(toFeat)
        if isfield(response.docs, toFeat{i})
            if ischar(response.docs.(toFeat{i}))
                tableTmp.(toFeat{i}) = {response.docs.(toFeat{i})};
            else
                nElem = length(response.docs.(toFeat{i}));
                if nElem > 1
                    strEle = '';
                    for k=1:nElem
                        if isnumeric(response.docs.(toFeat{i})(k))
                            strEle = strcat(strEle, num2str(response.docs.(toFeat{i})(k)), ';');
                        else
                            strEle = strcat(strEle, response.docs.(toFeat{i})(k), ';');
                        end
                    end
                    response.docs.(toFeat{i}) = {strEle};
                end
            end
            if iscell(response.docs.(toFeat{i}))
                tableTmp.(toFeat{i}) = response.docs.(toFeat{i});
            elseif isnumeric(response.docs.(toFeat{i}))
                tableTmp.(toFeat{i}) = response.docs.(toFeat{i});
            end
        else
            tableTmp.(toFeat{i}) = {'not found'};
        end
    end
end
end