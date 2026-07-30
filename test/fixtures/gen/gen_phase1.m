function gen_phase1(outdir)
% Generate golden-master fixtures for Phase 1 (Utils + Auxiliary) and the
% OLS engine that Phase 2 builds on.
%
% Written as a FUNCTION, not a script: the Replic/*/GO_*.m scripts open with
% `clear all`, which silently wipes a driver's loop state. Keeping generators
% function-scoped avoids that class of failure entirely.
%
% Usage:  gen_phase1('/abs/path/to/test/fixtures')
%
% VAR Toolbox 3.1 -> VARToolbox.jl port
% Fixture generator. See PROVENANCE.md for the reference environment.

if ~exist('outdir','var'); error('gen_phase1: outdir is required'); end
addpath(genpath('/Users/hiro/GitHub/VAR-Toolbox/v3dot1'));

F = struct();

%% ---- deterministic input blocks -------------------------------------
% Fixed, human-readable inputs. No RNG anywhere in this file.
A = [ 1.0  2.0  3.0
      4.0 -5.0  6.0
     -7.0  8.0  9.5
      2.5  0.5 -1.5
      3.3  4.4  5.5
     -0.2  1.1  2.2
      6.0 -6.0  0.0
      1.7  2.9 -3.1 ];
F.A = A;

% NaN patterns: leading run, trailing run, and an INTERIOR NaN. The interior
% one is what exposes CommonSample's `lo` behaviour (see PORTING_LOG).
N = [ NaN  NaN
      NaN  2.0
      1.0  3.0
      4.0  NaN     % <- interior NaN row
      5.0  6.0
      7.0  8.0
      NaN  NaN
      NaN  NaN ];
F.N = N;

[F.cs_out, F.cs_fo, F.cs_lo] = CommonSample(N);
[F.cs2_out, F.cs2_fo, F.cs2_lo] = CommonSample(N', 2);

% A clean matrix through CommonSample (no NaNs at all)
[F.cs3_out, F.cs3_fo, F.cs3_lo] = CommonSample(A);

%% ---- lag / design-matrix helpers ------------------------------------
F.mkl_2 = VARmakelags(A,2);
F.mkl_1 = VARmakelags(A,1);
for c = 0:3
    [Y,X] = VARmakexy(A,2,c);
    F.(sprintf('mkxy_Y_c%d',c)) = Y;
    F.(sprintf('mkxy_X_c%d',c)) = X;
end

F.lag_1   = lag(A,1);
F.lag_2   = lag(A,2);
F.lag_2v  = lag(A,2,-99);
F.trimr_a = trimr(A,2,1);
F.L_pos   = L(A,2);    % NB: positive nlag is a LEAD in L.m
F.L_neg   = L(A,-2);
F.L_zero  = L(A,0);
F.seqa    = seqa(2,0.5,7);
F.ptrend_0 = ptrend(0,8);
F.ptrend_2 = ptrend(2,8);
F.getqr    = getqr(A(1:3,1:3));
F.getqr_b  = getqr([2 -1 0; -1 2 -1; 0 -1 2]);

F.num2nan = Num2NaN([123456789 2; 3 123456789; 5 6]);
F.nan2num = NaN2Num([NaN 2; 3 NaN; 5 6]);

%% ---- prctile: the M1 check ------------------------------------------
% MATLAB prctile is the Hazen/type-5 definition. Julia's default quantile is
% type 7. These fixtures pin the difference.
pv = [3.1 1.4 7.7 2.2 9.6 5.5 4.4 8.8]';
F.pv = pv;
F.pct_levels = [0 2.5 25 50 75 95 97.5 100]';
F.pct_vec = prctile(pv, F.pct_levels);
F.pct_mat_d1 = prctile(A, [2.5 50 97.5], 1);
F.pct_mat_d2 = prctile(A, [2.5 50 97.5], 2);
% 4-D case, exercised by SR.m / VARirband.m with dim=4
R4 = reshape(linspace(-3, 4, 3*4*2*9), [3 4 2 9]);
F.R4 = R4;
F.pct_4d = prctile(R4, [16 50 84], 4);
% NaN handling: prctile omits NaNs
pn = pv; pn([2 5]) = NaN;
F.pn = pn;
F.pct_nan = prctile(pn, F.pct_levels);

%% ---- distribution functions: the M3 check ---------------------------
% These are hand-rolled on base MATLAB betainc, NOT the Statistics Toolbox.
% tdis_prb in particular clamps x2 to 1-1e-12, which perturbs the value at
% x=0 at the 1e-6 level because the Beta(n/2,1/2) density is singular at 1.
tx = [-4.5 -2.0 -1.0 -0.25 0 0.25 1.0 2.0 4.5 25.0]';
F.tx = tx;
F.tdis_prb_5   = tdis_prb(tx, 5);
F.tdis_prb_30  = tdis_prb(tx, 30);
F.tdis_prb_200 = tdis_prb(tx, 200);

tp = [0.001 0.025 0.1 0.25 0.5 0.75 0.9 0.975 0.999]';
F.tp = tp;
F.tdis_inv_5   = tdis_inv(tp, 5);
F.tdis_inv_30  = tdis_inv(tp, 30);
F.tdis_inv_200 = tdis_inv(tp, 200);

F.beta_pdf_23 = beta_pdf([0.01 0.1 0.5 0.9 0.99]', 2, 3);
F.beta_inv_23 = beta_inv([0.01 0.1 0.5 0.9 0.99]', 2, 3);
F.fdis_prb_35 = fdis_prb([0.1 0.5 1 2 5 10]', 3, 5);

%% ---- nan* family ----------------------------------------------------
% Ported to Statistics/StatsBase equivalents; these pin the conventions
% (notably the n-1 denominator in nanstd/nanvar and the all-NaN -> NaN rule).
M = [ 1.0  NaN  3.0
      4.0  5.0  NaN
      NaN  NaN  NaN
      2.0  8.0  6.0 ];
F.M = M;
F.nanmean_1 = nanmean(M,1);   F.nanmean_2 = nanmean(M,2);
F.nansum_1  = nansum(M,1);    F.nansum_2  = nansum(M,2);
F.nanstd_1  = nanstd(M,1);    F.nanvar_1  = nanvar(M,1);
F.nanmedian_1 = nanmedian(M,1);
F.nansem_1  = nansem(M,1);

%% ---- OLSmodel -------------------------------------------------------
y = A(:,1);
x = A(:,2:3);
for c = 0:2   % NB: const==3 is BROKEN upstream (see below) -- excluded
    O = OLSmodel(y,x,c);
    pre = sprintf('ols_c%d_',c);
    F.([pre 'beta'])   = O.beta;
    F.([pre 'yhat'])   = O.yhat;
    F.([pre 'resid'])  = O.resid;
    F.([pre 'sige'])   = O.sige;
    F.([pre 'bstd'])   = O.bstd;
    F.([pre 'tstat'])  = O.tstat;
    F.([pre 'tprob'])  = O.tprob;
    F.([pre 'bint'])   = O.bint;
    F.([pre 'rsqr'])   = O.rsqr;
    F.([pre 'rbar'])   = O.rbar;
    F.([pre 'dw'])     = O.dw;
    F.([pre 'nvar'])   = O.nvar;
    F.([pre 'nobs'])   = O.nobs;
    F.([pre 'sigbeta'])= O.sigbeta;
    F.([pre 'bstd_HW'])= O.bstd_HW;
    F.([pre 'bstd_NW'])= O.bstd_NW;
    if c > 0; F.([pre 'F']) = O.F; end
end

% Record that const==3 raises. OLSmodel.m adds only 2 columns
% ([ones trend.^2]) but increments nvar by 3, so xpxi is built with the
% wrong dimension. VARmakexy.m's const==3 adds 3 columns and is fine --
% the two disagree.
F.ols_c3_errors = 0;
try
    OLSmodel(y,x,3);
catch
    F.ols_c3_errors = 1;
end

%% ---- save + report --------------------------------------------------
if ~exist(outdir,'dir'); mkdir(outdir); end
out = fullfile(outdir,'phase1.mat');
save(out,'-struct','F','-v7');
fprintf('wrote %s (%d fields)\n', out, numel(fieldnames(F)));
end
