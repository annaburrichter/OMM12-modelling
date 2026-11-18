function run_rfastcormics_batch()
%% run_rfastcormics_batch
%
% Runs rFASTCORMICS v2 for 12 species × 4 conditions using metaproteomics
% gene lists. Produces one context-specific SBML model per combination.
%
% Call from MATLAB command window: run_rfastcormics_batch()
%
% Assumptions
%   - Gene IDs in the SBML models and in the *_genes.csv files share the
%     same format (e.g. QQQ91203.1) — no translation needed.
%   - Detected gene -> 1, all other genes -> 0 (no hard-absent, -1)
%   - The "mono" condition file is named <SPECIES>-<SPECIES>_genes.csv.
%   - SBML models are gzip-compressed XML: <SPECIES>.xml.gz.
%
% Requirements
%   - COBRA Toolbox (initCobraToolbox)
%   - rFASTCORMICS v2 (https://github.com/sysbiolux/rFASTCORMICS)
%     its scripts/ folder must be added to the path AFTER COBRA so that its
%     local fastcore shadows COBRA's version.
%   - A compatible LP solver (Gurobi or CPLEX).

% =========================================================================
%  USER SETTINGS  –  edit these before running
% =========================================================================

MODEL_DIR   = 'OMM_munich_matlab';          % folder containing <SPECIES>.xml.gz files
GENE_DIR    = 'OMM_munich_matlab';          % folder containing *_genes.csv files
OUTPUT_DIR  = 'OMM_munich_matlab/csm_rfastcormics';  % where output models are written

RFASTCORMICS_SCRIPTS = 'software/rFASTCORMICS/rFASTCORMICS for RNA-seq data/rFASTCORMICS_v2/scripts';
SOLVER      = 'gurobi';     % 'gurobi' or 'ibm_cplex'

% rFASTCORMICS parameters
EPSILON              = 1e-4;     % FASTCORE flux threshold
FILLING_MEDIUM_FLAG  = 1;        % 1 = try to fill gaps with medium reactions
ADAPTIVE_SCALING     = 0;        % 0 recommended for proteomics data

% Biomass reaction name. If left empty ('') the script auto-detects it per
% model by searching for 'biomass' or 'growth' (case-insensitive) in rxns.
BIOMASS_NAME = 'bio1';   % e.g. 'BIOMASS_reaction'

% Reactions that must appear in every context-specific model (beyond biomass).
FUNC_RXNS = {};

% Exchange metabolite IDs that define the medium. Leave empty {} to use the
% bounds already encoded in the model.
MEDIUM_METS = {};

% Per-species biomass override (takes priority over BIOMASS_NAME)
BIOMASS_OVERRIDE = containers.Map({'YL32'}, {'bio1'});

% YL32 needs a smaller epsilon
EPSILON_OVERRIDE = containers.Map({'YL32'}, {1e-6});

% Run only these species ({} = all). Useful for reruns.
RUN_ONLY = {''};

% true = rebuild models even if output files already exist
OVERWRITE = true;

% =========================================================================
%  SPECIES & CONDITIONS
% =========================================================================

SPECIES = {'I46','I48','I49','KB18','KB1','YL27','YL2','YL31','YL32','YL44','YL45','YL58'};
CONDITIONS = {'AAM','M10','mouse','mono'};
if ~isempty(RUN_ONLY)
    SPECIES = SPECIES(ismember(SPECIES, RUN_ONLY));
end

% =========================================================================
%  SETUP
% =========================================================================

if ~exist(OUTPUT_DIR, 'dir'), mkdir(OUTPUT_DIR); end

logfile = fullfile(OUTPUT_DIR, 'batch_log.txt');
flog    = fopen(logfile, 'a');

logmsg(flog, 'Batch started.');

fprintf('Initialising COBRA Toolbox...\n');
initCobraToolbox(false);
changeCobraSolver(SOLVER, 'LP');

% Add rFASTCORMICS AFTER COBRA so its fastcore shadows COBRA's version.
addpath(genpath(RFASTCORMICS_SCRIPTS));
verifyFastcorePath();

feature astheightlimit 2000;   % raise MATLAB stack limit for recursive calls

% =========================================================================
%  MAIN LOOP
% =========================================================================

totalModels  = numel(SPECIES) * numel(CONDITIONS);
doneCount    = 0;
skippedCount = 0;
failedList   = {};

for si = 1:numel(SPECIES)
    sp = SPECIES{si};

    %% Load generic model once per species ----------------------------
    gz_file  = fullfile(MODEL_DIR, [sp '.xml.gz']);

    if ~exist(gz_file, 'file')
        warning('Model file not found: %s — skipping species %s.', gz_file, sp);
        logmsg(flog, 'SKIP species %s — model file not found.', sp);
        skippedCount = skippedCount + numel(CONDITIONS);
        continue
    end

    fprintf('\n======= %s =======\n', sp);
    logmsg(flog, 'Loading model for species %s', sp);

    % Decompress to a temp location, load, then delete
    gunzip(gz_file, OUTPUT_DIR);
    unzipped = fullfile(OUTPUT_DIR, [sp '.xml']);
    model_base = readCbModel(unzipped);
    delete(unzipped);

    % rFASTCORMICS v2 does vertcat(model.subSystems{:}), so every entry must be
    % a 1x1 cell containing a char row vector, e.g. {'Glycolysis'}.
    if ~isfield(model_base, 'subSystems')
        model_base.subSystems = repmat({''}, numel(model_base.rxns), 1);
    end
    ss = model_base.subSystems;
    for k = 1:numel(ss)
        s = ss{k};
        if iscell(s)                       % {} / {'A'} / {'A','B'}
            s = strjoin(cellfun(@char, s(:)', 'UniformOutput', false), ';');
        end
        ss{k} = {char(s)};                 % '' -> {''}, 'A' -> {'A'}
    end
    model_base.subSystems = ss(:);

    % Resolve biomass reaction name for this species
    if isKey(BIOMASS_OVERRIDE, sp)
        bmSetting = BIOMASS_OVERRIDE(sp);
    else
        bmSetting = BIOMASS_NAME;
    end
    bm = resolveBiomassName(model_base, bmSetting, sp);
    if isempty(bm)
        warning('Cannot determine biomass reaction for %s — skipping.', sp);
        logmsg(flog, 'SKIP species %s — biomass reaction not found.', sp);
        skippedCount = skippedCount + numel(CONDITIONS);
        continue
    end
    model_base = changeObjective(model_base, bm);
    fprintf('  Biomass: %s  |  Genes in model: %d  |  Rxns: %d\n', ...
            bm, numel(model_base.genes), numel(model_base.rxns));

    %% Loop over conditions -------------------------------------------
    for ci = 1:numel(CONDITIONS)
        cond = CONDITIONS{ci};

        % mono gene file is named <SPECIES>-<SPECIES>_genes.csv
        if strcmp(cond, 'mono')
            gene_suffix = sp;
        else
            gene_suffix = cond;
        end
        gene_file = fullfile(GENE_DIR, [sp '-' gene_suffix '_genes.csv']);

        out_sbml  = fullfile(OUTPUT_DIR, [sp '_' cond '.xml']);
        out_stats = fullfile(OUTPUT_DIR, [sp '_' cond '_stats.txt']);

        fprintf('  [%s] ', cond);

        % Skip if output already exists
        if exist(out_sbml, 'file') && ~OVERWRITE
            fprintf('already exists — skipping.\n');
            logmsg(flog, '%s | %s — output exists, skipped.', sp, cond);
            skippedCount = skippedCount + 1;
            continue
        end

        % Check gene file
        if ~exist(gene_file, 'file')
            warning('Gene file not found: %s — skipping.', gene_file);
            logmsg(flog, '%s | %s — gene file not found: %s', sp, cond, gene_file);
            skippedCount = skippedCount + 1;
            continue
        end

        %% Build discretized vector & dico --------------------------------
        try
            detectedGenes = readDetectedGenes(gene_file);
        catch ME
            warning('Failed to read %s: %s', gene_file, ME.message);
            logmsg(flog, '%s | %s — gene file read error: %s', sp, cond, ME.message);
            failedList{end+1} = sprintf('%s|%s', sp, cond); %#ok<AGROW>
            continue
        end

        genes  = model_base.genes;   % IDs as stored by COBRA (may be mangled)
        nGenes = numel(genes);

        % Normalize model gene IDs to match CSV format.
        % COBRA's SBML loader prepends 'G_' and replaces '.' and '-' with '_'.
        % We reverse this so the IDs align with the metaproteomics CSV.
        % normalizeModelGenes() tries several transformations and picks the
        % one with the highest overlap against detectedGenes.
        [normGenes, transform, overlap] = normalizeModelGenes(genes, detectedGenes);

        if overlap == 0
            % Print diagnostics to help the user fix manually
            fprintf('\n  [DIAG] 0 gene overlap for %s|%s.\n', sp, cond);
            fprintf('  Model genes (first 5):    %s\n', strjoin(genes(1:min(5,end)),  ', '));
            fprintf('  CSV genes   (first 5):    %s\n', strjoin(detectedGenes(1:min(5,end)), ', '));
            fprintf('  Normalized  (first 5):    %s\n', strjoin(normGenes(1:min(5,end)), ', '));
            warning('No detected genes match model genes for %s|%s — skipping.', sp, cond);
            logmsg(flog, '%s | %s — 0 gene overlap (transform: %s), skipped.', sp, cond, transform);
            skippedCount = skippedCount + 1;
            continue
        end
        fprintf('genes in model: %d / %d detected (transform: %s)   ', ...
                overlap, numel(detectedGenes), transform);

        discretized = zeros(nGenes, 1);           % 0 = unknown
        discretized(ismember(normGenes, detectedGenes)) = 1;

        % dico maps model gene IDs (COBRA format) → data gene IDs (CSV format)
        dico = table(genes, normGenes, 'VariableNames', {'modelGene', 'dataGene'});

        %% Build optionalSettings -----------------------------------------
        opt = struct();
        func_rxns = FUNC_RXNS;
        if ~ismember(bm, func_rxns)
            func_rxns{end+1} = bm; %#ok<AGROW>
        end
        opt.func = func_rxns;
        if ~isempty(MEDIUM_METS)
            opt.medium = MEDIUM_METS;
        end

        %% Run rFASTCORMICS -----------------------------------------------
        try
            eps_sp = EPSILON;
            if isKey(EPSILON_OVERRIDE, sp), eps_sp = EPSILON_OVERRIDE(sp); end
            [csm, ~, coreIdx] = rFastcormics( ...
                model_base, discretized, genes, dico, ...
                bm, 1, eps_sp, opt, ...
                FILLING_MEDIUM_FLAG, ADAPTIVE_SCALING);
        catch ME
            fprintf(2, '\n%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
            warning('rFastcormics failed for %s|%s: %s', sp, cond, ME.message);
            logmsg(flog, '%s | %s — rFastcormics error: %s', sp, cond, ME.message);
            failedList{end+1} = sprintf('%s|%s', sp, cond); %#ok<AGROW>
            continue
        end

        %% Save SBML -------------------------------------------------------
        try
            csm = changeObjective(csm, bm);
            writeCbModel(csm, 'format', 'sbml', 'fileName', out_sbml);
        catch ME
            warning('writeCbModel failed for %s|%s: %s', sp, cond, ME.message);
            logmsg(flog, '%s | %s — writeCbModel error: %s', sp, cond, ME.message);
            failedList{end+1} = sprintf('%s|%s', sp, cond); %#ok<AGROW>
            continue
        end

        %% Write stats file ------------------------------------------------
        writeStats(out_stats, sp, cond, gene_file, overlap, numel(detectedGenes), ...
                   nGenes, numel(coreIdx), numel(csm.rxns), bm);

        doneCount = doneCount + 1;
        fprintf('-> CSM rxns: %d (core: %d)  OK\n', numel(csm.rxns), numel(coreIdx));
        logmsg(flog, '%s | %s — done. CSM rxns: %d, core: %d', ...
               sp, cond, numel(csm.rxns), numel(coreIdx));
    end
end

%% =========================================================================
%  SUMMARY
% =========================================================================

fprintf('\n========================================\n');
fprintf('Batch complete.\n');
fprintf('  Done:    %d / %d\n', doneCount, totalModels);
fprintf('  Skipped: %d\n',      skippedCount);
fprintf('  Failed:  %d\n',      numel(failedList));
if ~isempty(failedList)
    fprintf('  Failed combos:\n');
    for k = 1:numel(failedList)
        fprintf('    %s\n', failedList{k});
    end
end
fprintf('Models saved to: %s\n', OUTPUT_DIR);
logmsg(flog, 'Batch finished. Done=%d Skipped=%d Failed=%d', doneCount, skippedCount, numel(failedList));
fclose(flog);

end % main function


% =========================================================================
%  LOCAL HELPER FUNCTIONS
% =========================================================================

function logmsg(flog, fmt, varargin)
    fprintf(flog, ['[%s] ' fmt '\n'], datestr(now,'yyyy-mm-dd HH:MM:SS'), varargin{:});
end

function verifyFastcorePath()
% Warns if COBRA's fastcore is shadowing the rFASTCORMICS version.
    hits = which('fastcore', '-all');
    if isempty(hits)
        error('fastcore not found. Add rFASTCORMICS scripts to the path.');
    end
    if ~contains(hits{1}, 'rFASTCORMICS')
        warning(['COBRA fastcore may shadow rFASTCORMICS fastcore.\n' ...
                 'Ensure addpath(rFASTCORMICS) comes AFTER initCobraToolbox.\n' ...
                 'First hit: %s'], hits{1});
    end
end

function bm = resolveBiomassName(model, bmSetting, spName)
% Returns the biomass reaction name for this species.
% Priority:
%   1. BIOMASS_NAME setting (if set and found in model.rxns)
%   2. Name-pattern match: biomass|growth|objective|maintenance
%   3. Objective coefficient vector (model.c ~= 0)  <-- most reliable fallback
% Prints candidate names to help diagnose misses.

    % --- 1. Explicit setting ---
    if ~isempty(bmSetting)
        if ismember(bmSetting, model.rxns)
            bm = bmSetting;
            return
        else
            warning('BIOMASS_NAME "%s" not in model for %s. Auto-detecting.', bmSetting, spName);
        end
    end

    % --- 2. Name-pattern match ---
    idx = find(~cellfun(@isempty, regexpi(model.rxns, ...
        '(biomass|growth|objective|maintenance|atpm)', 'match')));
    if numel(idx) == 1
        bm = model.rxns{idx};
        fprintf('  [auto-detected biomass by name: %s]\n', bm);
        return
    elseif numel(idx) > 1
        % Prefer one that contains 'biomass'
        bm_sub = find(~cellfun(@isempty, regexpi(model.rxns(idx), 'biomass', 'match')));
        if numel(bm_sub) >= 1
            bm = model.rxns{idx(bm_sub(1))};
        else
            bm = model.rxns{idx(1)};
        end
        fprintf('  [auto-detected biomass by name: %s (from %d candidates: %s)]\n', ...
                bm, numel(idx), strjoin(model.rxns(idx), ', '));
        return
    end

    % --- 3. Objective coefficient vector ---
    if isfield(model, 'c') && any(model.c ~= 0)
        obj_idx = find(model.c ~= 0);
        if numel(obj_idx) == 1
            bm = model.rxns{obj_idx};
            fprintf('  [auto-detected biomass via objective vector: %s]\n', bm);
        else
            bm = model.rxns{obj_idx(1)};
            fprintf('  [auto-detected biomass via objective vector (first of %d): %s]\n', ...
                    numel(obj_idx), bm);
        end
        return
    end

    % --- Nothing found: print all rxns to help diagnose ---
    bm = '';
    fprintf('  [ERROR] Could not detect biomass for %s.\n', spName);
    fprintf('  First 20 reaction IDs in model:\n');
    for k = 1:min(20, numel(model.rxns))
        fprintf('    %s\n', model.rxns{k});
    end
    fprintf('  Set BIOMASS_NAME explicitly in the USER SETTINGS block.\n');
end

function [normGenes, bestTransform, bestOverlap] = normalizeModelGenes(modelGenes, csvGenes)
% Tries several transformations to reverse COBRA's SBML gene ID mangling,
% picks the one that gives the highest overlap with csvGenes.
%
% COBRA SBML loading converts e.g. QQQ91203.1 -> G_QQQ91203_1
% (prepends 'G_', replaces '.' and '-' with '_').
%
% Transformations tried (in order):
%   'none'        - as-is
%   'stripG'      - remove leading 'G_'
%   'dot'         - remove 'G_', replace last _N suffix with .N
%   'dotAll'      - remove 'G_', replace ALL underscores with dots
%   'dotDash'     - remove 'G_', replace underscores with dots except
%                   internal ones that were dashes

    transforms = {'none', 'stripG', 'dot', 'dotAll'};
    bestOverlap   = -1;
    bestTransform = 'none';
    normGenes     = modelGenes;

    for t = 1:numel(transforms)
        tr = transforms{t};
        cand = applyTransform(modelGenes, tr);
        ov = sum(ismember(cand, csvGenes));
        if ov > bestOverlap
            bestOverlap   = ov;
            bestTransform = tr;
            normGenes     = cand;
        end
        if ov == numel(csvGenes), break; end   % perfect match, stop early
    end
end

function out = applyTransform(genes, tr)
    switch tr
        case 'none'
            out = genes;
        case 'stripG'
            % Remove leading 'G_'
            out = regexprep(genes, '^G_', '');
        case 'dot'
            % Remove 'G_', then replace the last _<digits> with .<digits>
            % e.g. G_QQQ91203_1 -> QQQ91203.1
            out = regexprep(genes, '^G_', '');
            out = regexprep(out,   '_(\d+)$', '.$1');
        case 'dotAll'
            % Remove 'G_', replace all underscores with dots
            out = regexprep(genes, '^G_', '');
            out = strrep(out, '_', '.');
        otherwise
            out = genes;
    end
end

function genes = readDetectedGenes(csvfile)
% Reads a single-column CSV (header containing 'gene') and returns a cell array.
    T = readtable(csvfile, 'TextType', 'string');
    colnames = lower(T.Properties.VariableNames);
    col = find(contains(colnames, 'gene'), 1);
    if isempty(col)
        error('No column containing "gene" found in %s. Columns: %s', ...
              csvfile, strjoin(T.Properties.VariableNames, ', '));
    end
    genes = cellstr(T{:, col});
    genes = strtrim(genes);
    genes = genes(~cellfun(@isempty, genes));
end

function writeStats(outfile, sp, cond, gene_file, overlap, nDetected, nModelGenes, coreSize, modelSize, bm)
% Writes a short text summary for one species-condition pair.
    fid = fopen(outfile, 'w');
    fprintf(fid, 'Species:             %s\n', sp);
    fprintf(fid, 'Condition:           %s\n', cond);
    fprintf(fid, 'Gene file:           %s\n', gene_file);
    fprintf(fid, 'Biomass reaction:    %s\n', bm);
    fprintf(fid, 'Detected proteins:   %d\n', nDetected);
    fprintf(fid, 'Genes in model:      %d\n', nModelGenes);
    fprintf(fid, 'Overlap:             %d (%.1f%% of detected in model)\n', ...
            overlap, 100 * overlap / nDetected);
    fprintf(fid, 'Core reactions:      %d\n', coreSize);
    fprintf(fid, 'CSM reactions:       %d\n', modelSize);
    fprintf(fid, 'Generated:           %s\n', datestr(now));
    fclose(fid);
end
