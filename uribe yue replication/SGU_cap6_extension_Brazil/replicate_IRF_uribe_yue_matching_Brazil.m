%% Uribe-Yue: IRF matching -- DSGE vs. VAR, and estimated vs. paper parameters
% ****************************************************************
% *** BRAZIL ONLY. Sample: 1996:Q1 - 2026:Q1 (single-country VAR ***
% *** replication, uribe_yue_VAR_Brazil.mlx / Table 6.1.BR).     ***
% *** NOT the paper's original 7-country panel (1994-2001) --     ***
% *** every result plotted and printed below reflects Brazil's    ***
% *** own, much longer and more current, estimated dynamics.      ***
% ****************************************************************
%
% Runs the IRF-matching Dynare model (DSGE_uribe_yue_IRFmatch_Brazil.mod),
% then plots the model-implied IRFs (post-estimation, from oo_.irfs)
% against the Brazil-only VAR target IRFs (from
% IRF_target_uribe_yue_Brazil.mat) as two lines per chart, and prints the
% estimated {PSSI, PHI, ETA, MU} against the values reported in the paper
% (Uribe and Yue, 2006, JIE, Table 3, 7-country panel 1994-2001) -- NOT
% against a Brazil-specific benchmark, since the paper never estimated one;
% the comparison here is "how do Brazil-implied frictions differ from the
% original panel's," not a replication check.
%
% Matches the paper's own IRF-matching specification (Section 3.1.4): 4
% variables (output, investment, trade balance/GDP, country rate), 2
% shocks, 24 quarters, both shocks sized at 0.01 ("unit innovation" / "one
% percentage point increase", per Figures 2-3's captions) -- the ONLY thing
% that changed from the baseline run is which country's estimated R and
% Rus processes (eqs. 6.20/6.22) and target IRFs feed the estimation. The
% world interest rate, rus, is NOT one of the paper's matched variables,
% so it is not plotted here -- only the 4 matched variables the paper
% itself compares.
%
% PREREQUISITE: run uribe_yue_VAR_Brazil.mlx through Section 3, then paste
% in and run uribe_yue_VAR_Brazil_target_export.m (Section 3c) to produce
% IRF_target_uribe_yue_Brazil.mat, BEFORE running this script.

%% Step 1: Run the IRF-matching Dynare model (Brazil-only law of motion)

clear
close all
clc

dynare DSGE_uribe_yue_IRFmatch_Brazil

%% Step 2: Load the Brazil-only VAR target IRFs
% (Also loaded internally by the .mod file itself, but re-loaded here
% explicitly so this script is self-contained even if run in a fresh
% session against saved results.)

if ~isfile('IRF_target_uribe_yue_Brazil.mat')
    error(['IRF_target_uribe_yue_Brazil.mat not found. Run uribe_yue_VAR_Brazil.mlx ' ...
           'through Section 3, then uribe_yue_VAR_Brazil_target_export.m (Section 3c), first.']);
end
target = load('IRF_target_uribe_yue_Brazil.mat');

HORIZON = 24;   % matched horizon, per the paper ("the first 24 quarters")
if isfield(target,'HORIZON'), HORIZON = target.HORIZON; end

%% Step 3: Plot -- DSGE (model, post-estimation) vs. VAR (target), two lines per
% chart, 2x2 grid: yy, tby, r, iv -- the paper's 4 matched variables only.

matched_vars = {'yy','tby','r','iv'};
shock_names  = {'eps_r','eps_rus'};
shock_titles = {'eps^r (country spread, unit/0.01 shock)','eps^{rus} (world rate, unit/0.01 shock)'};

for is_ = 1:numel(shock_names)
    sname = shock_names{is_};
    figure('Color','w','Name',['BRAZIL 1996-2026: DSGE vs VAR -- ' shock_titles{is_}]);

    for iv_ = 1:numel(matched_vars)
        vname = matched_vars{iv_};
        dsge_field   = [vname '_' sname];
        target_field = [vname '_' sname];

        subplot(2, 2, iv_)
        hold on

        if isfield(oo_.irfs, dsge_field)
            dsge_series = oo_.irfs.(dsge_field)(:);
            plot(0:numel(dsge_series)-1, dsge_series, 'b-', 'LineWidth', 1.75);
        else
            warning(['oo_.irfs.' dsge_field ' not found -- skipping DSGE line for ' vname]);
            dsge_series = [];
        end

        if isfield(target, target_field)
            target_series = target.(target_field)(:);
            plot(0:numel(target_series)-1, target_series, 'r--', 'LineWidth', 1.75);
        else
            warning(['target.' target_field ' not found -- skipping VAR line for ' vname]);
        end
        legend({'DSGE (estimated)','VAR target'}, 'Location', 'best')

        plot([0 numel(dsge_series)-1], [0 0], 'k-', 'LineWidth', 0.5)
        title(vname, 'Interpreter', 'none')
        xlabel('quarters after shock')
        ylabel('deviation from steady state')
        axis tight
        grid on
        hold off
    end
    annotation('textbox', [0 0.94 1 0.05], 'String', ...
        ['BRAZIL, 1996:Q1-2026:Q1 -- IRF matching: DSGE vs VAR -- ' shock_titles{is_}], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

%% Step 4: Print estimated parameters, compared to the paper

param_names_to_check = {'PSSI','PHI','ETA','MU'};
% Table 3 of Uribe and Yue (2006, JIE): psi=0.00042, phi=72.8, eta=1.2, mu=0.204.
% uribe_yue_run.m's own calibration (0.000425, 72.8268, 1.2023, 0.2037) carries
% more decimal places and is used here as the comparison benchmark. THIS IS
% STILL THE ORIGINAL 7-COUNTRY-PANEL (1994-2001) BENCHMARK -- the paper
% never estimated a Brazil-specific version of Table 3, so what follows is
% "how do Brazil's implied frictions differ from the original panel's,"
% not a check against a Brazil-specific published answer.
paper_values = [0.000425, 72.8268, 1.2023, 0.2037];

fprintf('\n============================================================\n');
fprintf(' BRAZIL ONLY, 1996:Q1-2026:Q1\n');
fprintf(' IRF-matching estimation results vs. the paper''s ORIGINAL panel\n');
fprintf(' (Uribe and Yue, 2006, JIE, Table 3, 7-country panel, 1994-2001)\n');
fprintf('============================================================\n');
fprintf('%-8s %14s %14s %12s\n', 'Param', 'Estimated (BR)', 'Paper (panel)', '% diff');

for j = 1:numel(param_names_to_check)
    pname = param_names_to_check{j};
    idx = find(strcmp(cellstr(M_.param_names), pname));
    if isempty(idx)
        warning(['Could not find parameter ' pname ' in M_.param_names.']);
        continue
    end
    est_value = M_.params(idx);
    paper_value = paper_values(j);
    pct_diff = 100*(est_value - paper_value)/paper_value;
    fprintf('%-8s %14.6f %14.6f %11.1f%%\n', pname, est_value, paper_value, pct_diff);
end
fprintf('============================================================\n');
fprintf(['NOTE: this estimation uses EQUAL weighting across all 192 target\n' ...
         'points, while the paper (eq. 17) inverse-variance weights by each\n' ...
         'point''s delta-method standard error. Differences from Table 3 above\n' ...
         'may partly reflect that, or the very different sample (single\n' ...
         'country, 1996-2026) and estimated R/Rus processes, rather than any\n' ...
         'error in the model itself.\n']);

% If available, also report method_of_moments' own diagnostic structure
% (field names can vary slightly across Dynare versions, so this is
% wrapped defensively and simply skipped if not present).
try
    if isfield(oo_, 'mom')
        fprintf('\noo_.mom is also available for further diagnostics (e.g. standard errors,\n');
        fprintf('objective function value at the optimum) -- inspect oo_.mom directly.\n');
    end
catch
    % no-op: oo_.mom structure not present in this Dynare version/run
end

%% Step 5: Save everything for later reference

save('uribe_yue_IRFmatch_results_Brazil.mat', 'oo_', 'M_', 'target', 'param_names_to_check', 'paper_values');
fprintf('\nSaved estimation results to uribe_yue_IRFmatch_results_Brazil.mat\n');
