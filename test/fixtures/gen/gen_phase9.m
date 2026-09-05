function gen_phase9(outdir)
% Golden-master fixtures for Phase 9: upstream v4.0 bug fixes adopted by
% the port.
%
% Unlike gen_phase1..8, which run VAR Toolbox 3.1 (v3dot1/), this generator
% runs VAR Toolbox 4.0 (v4dot0/ -- upstream ambropo/VAR-Toolbox@36e706a,
% "v4.0 release" + "Fixed a few minor bugs"). v4.0 redesigned the API
% (VARmodel now takes VARopt as its 4th input and performs identification
% itself; VARir/VARvd/VARhd/VAR*band no longer exist), so each block below
% drives the v4.0 entry point and extracts the quantity whose v3.1
% computation was defective and which the port now computes the v4.0 way:
%
%   * OLSmodel / VARmodel / ARDLmodel  -- residual degrees of freedom in
%     tprob and bint (v3.1: nobs, or nobs - lags only)
%   * CommonSample                     -- lo is the TRAILING NaN count
%   * compute_IR                       -- shut forces the companion recursion
%   * compute_HD                       -- lagged exogenous regressors enter
%   * recover_B (iv)                   -- date-mask sample alignment
%   * LPmodel                          -- intercept present with no controls
%   * DatesPlot (monthly)              -- exact 1/12 arithmetic
%
% Inputs are read back from the v3.1 fixtures (phase1..5.mat) so that the
% two generations share bit-identical data.
%
% Usage:  gen_phase9('/abs/path/to/test/fixtures')

if ~exist('outdir','var'); error('gen_phase9: outdir is required'); end
root = '/Users/hiro/GitHub/VAR-Toolbox/v4dot0';
addpath(genpath(root));
cleanup = onCleanup(@() rmpath(genpath(root)));

F = struct();
F.source = 'v4dot0 (ambropo/VAR-Toolbox@36e706a)';

P1 = load(fullfile(outdir,'phase1.mat'), 'A', 'N');
P2 = load(fullfile(outdir,'phase2.mat'), 'A', 'E');
P3 = load(fullfile(outdir,'phase3.mat'), 'SW');
P4 = load(fullfile(outdir,'phase4.mat'), 'gk_ENDO', 'gk_IV');
P5 = load(fullfile(outdir,'phase5.mat'), 'SW');

%% ---- OLSmodel: residual degrees of freedom -----------------------------
y = P1.A(:,1); x = P1.A(:,2:3);
for c = 0:2
    O = OLSmodel(y,x,c);
    p = sprintf('ols_c%d_', c);
    F.([p 'tprob']) = O.tprob;
    F.([p 'bint'])  = O.bint;
    F.([p 'nobs'])  = O.nobs;
    F.([p 'nvar'])  = O.nvar;
end

%% ---- CommonSample: trailing count, empty input ------------------------
N = P1.N;
[OUT, fo, lo] = CommonSample(N);
F.cs_out = OUT; F.cs_fo = fo; F.cs_lo = lo;
[OUT2, fo2, lo2] = CommonSample(N', 2);
F.cs2_out = OUT2; F.cs2_fo = fo2; F.cs2_lo = lo2;
[OUTe, foe, loe] = CommonSample([]);
F.cs_empty_ok = double(isempty(OUTe) && foe == 0 && loe == 0);

%% ---- VARmodel: tprob with residual degrees of freedom -----------------
o = VARoption; o.ident = '';           % reduced form only
A = P2.A; E = P2.E;
cfg = {2,1; 1,1; 3,1; 2,0; 2,2};       % gen_phase2's list minus const=3,
                                       % which v4.0 rejects
for k = 1:size(cfg,1)
    V = VARmodel(A, cfg{k,1}, cfg{k,2}, o);
    p = sprintf('v%d_', k);
    F.([p 'ntotcoeff']) = V.ntotcoeff;
    for j = 1:V.nvar
        F.([p sprintf('eq%d_tprob', j)]) = V.(V.eqnames{j}).tprob;
    end
end
ex = {2,2; 3,1; 1,3};
for k = 1:size(ex,1)
    V = VARmodel(A, ex{k,1}, 1, o, E, ex{k,2});
    p = sprintf('vx%d_', k);
    F.([p 'ntotcoeff']) = V.ntotcoeff;
    for j = 1:V.nvar
        F.([p sprintf('eq%d_tprob', j)]) = V.(V.eqnames{j}).tprob;
    end
end
SW = P3.SW;
V = VARmodel(SW, 4, 1, o);
for j = 1:V.nvar
    F.(sprintf('sw_eq%d_tprob', j)) = V.(V.eqnames{j}).tprob;
end

%% ---- ARDLmodel: residual degrees of freedom ---------------------------
y1 = A(:,1);
D = ARDLmodel(y1,2,1);
F.ardl_tprob = D.tprob; F.ardl_bint = D.bint;
D2 = ARDLmodel(y1,2,1,E,1);
F.ardl2_tprob = D2.tprob; F.ardl2_bint = D2.bint;

%% ---- compute_IR: shut under the default recurs = 'wold' ---------------
oo = o; oo.ident = 'short'; oo.shut = 2; oo.nsteps = 12; oo.inference = 0;
V = VARmodel(SW, 4, 1, oo);
F.ir_shut2 = V.IR;
oo.recurs = 'comp';
V = VARmodel(SW, 4, 1, oo);
F.ir_shut2_comp = V.IR;

%% ---- compute_HD: lagged exogenous regressors --------------------------
% Two deterministic pseudo-random exogenous series (same LCG as gen_phase2,
% different seed), as long as the Stock-Watson sample.
T = size(SW,1);
lcg = zeros(2*T,1); s = 24680;
for i = 1:2*T
    s = mod(1103515245*s + 12345, 2^31);
    lcg(i) = s / 2^31 - 0.5;
end
EXO = [lcg(1:T), lcg(T+1:2*T)];
F.hd_EXOG = EXO;
hdcfg = {4,1,1; 4,2,1; 4,2,2; 2,3,1};  % nlag, nlag_ex, const
for k = 1:size(hdcfg,1)
    oo = o; oo.ident = 'short'; oo.inference = 0; oo.nsteps = 12;
    V = VARmodel(SW, hdcfg{k,1}, hdcfg{k,3}, oo, EXO, hdcfg{k,2});
    p = sprintf('hd%d_', k);
    F.([p 'nlag']) = hdcfg{k,1}; F.([p 'nlag_ex']) = hdcfg{k,2};
    F.([p 'const']) = hdcfg{k,3};
    F.([p 'shock']) = V.HD.shock; F.([p 'init']) = V.HD.init;
    F.([p 'const_c']) = V.HD.const; F.([p 'trend']) = V.HD.trend;
    F.([p 'exo']) = V.HD.exo;     F.([p 'endo']) = V.HD.endo;
    % Reconstruction identity, recorded so the Julia test can assert the
    % same bound MATLAB achieves.
    F.([p 'recon_err']) = max(abs(V.HD.endo(V.nlag+1:end,:) - V.Y), [], 'all');
end

%% ---- recover_B (iv): instrument with trailing / interior gaps ---------
ENDO = P4.gk_ENDO; IV = P4.gk_IV;
ivcase = {'trail', 'gap', 'both'};
for k = 1:3
    IVm = IV;
    if any(strcmp(ivcase{k}, {'trail','both'})); IVm(end-23:end) = NaN; end
    if any(strcmp(ivcase{k}, {'gap','both'}));   IVm(250:255)    = NaN; end
    oo = o; oo.ident = 'iv'; oo.IV = IVm; oo.nsteps = 48; oo.inference = 0;
    V = VARmodel(ENDO, 12, 1, oo);
    p = ['iv_' ivcase{k} '_'];
    F.([p 'IV'])      = IVm;
    F.([p 'B1'])      = V.B(:,1);
    F.([p 'IR1'])     = V.IR(:,:,1);
    F.([p 'sigma_b']) = V.sigma_b;
    F.([p 'fsbeta'])  = V.FirstStage.beta;
    % v4.0 multiplies Biv by sign(FirstStage.beta(2)) -- an undocumented
    % convention change relative to 3.1 that the port does not adopt.
    % Recorded so the test can divide it out.
    F.([p 'fssign'])  = sign(V.FirstStage.beta(2));
    F.([p 'nkeep'])   = size(V.FirstStage.x, 1);
end

%% ---- LPmodel: intercept with no controls ------------------------------
LPopt = LPoption; LPopt.nsteps = 8; LPopt.impact = 0;
ENDO5 = P5.SW(:,1); S5 = P5.SW(:,3);
for c = 1:2
    LP = LPmodel(ENDO5, S5, [], 4, c, LPopt);
    p = sprintf('lp_c%d_', c);
    F.([p 'IR'])      = LP.IR;
    F.([p 'h1_beta']) = LP.h1.beta;
    F.([p 'h1_nobs']) = LP.h1.nobs;
    F.([p 'h1_x'])    = LP.h1.x;
    F.([p 'ntotcoeff']) = LP.ntotcoeff;
end

%% ---- DatesPlot: monthly ticks with exact 1/12 arithmetic --------------
dp = {1965, 100, 8; 1965+5/12, 131, 8; 1990+11/12, 50, 5; 1959+3/12, 396, 8};
for k = 1:size(dp,1)
    fig = figure('Visible','off');
    plot(1:dp{k,2});
    DatesPlot(dp{k,1}, dp{k,2}, dp{k,3}, 'm');
    xt = get(gca,'XTick');
    xl = cellfun(@str2double, cellstr(get(gca,'XTickLabel')));
    xlim_ = get(gca,'XLim');
    close(fig);
    p = sprintf('dp%d_', k);
    F.([p 'fo']) = dp{k,1}; F.([p 'nobs']) = dp{k,2}; F.([p 'nticks']) = dp{k,3};
    F.([p 'xt']) = xt(:); F.([p 'xl']) = xl(:); F.([p 'xlim']) = xlim_;
end

%% ---- Stats guards: PairCorr with non-positive levels ------------------
% Deterministic strictly positive levels (LCG again), then a negative and a
% zero entry. v4.0 warns and skips the log first differences; the resulting
% first-difference column is whatever (sum(NaN,'omitnan')-1)./(n-1) gives.
lcg = zeros(90,1); s = 13579;
for i = 1:90
    s = mod(1103515245*s + 12345, 2^31);
    lcg(i) = s / 2^31;                  % ~U(0,1)
end
XP = 1 + reshape(lcg, 30, 3);
F.pc_X = XP;
F.pc_pos = PairCorr(XP);
XN = XP; XN(7,2) = -0.5;
F.pc_Xneg = XN;
warning('off','all'); F.pc_neg = PairCorr(XN); warning('on','all');
XZ = XP; XZ(3,1) = 0;
F.pc_Xzero = XZ;
warning('off','all'); F.pc_zero = PairCorr(XZ); warning('on','all');
% MovAvgCent / MovCorrCent short-series guards (contract, not a number).
F.movavgcent_short_errors = 0; F.movcorrcent_short_errors = 0;
try, MovAvgCent(XP(1:10,:), 5);     catch, F.movavgcent_short_errors  = 1; end
try, MovCorrCent(XP(1:10,:), 5, 1); catch, F.movcorrcent_short_errors = 1; end
F.movavgcent_edge_ok = double(all(size(MovAvgCent(XP(1:9,:), 4)) == [9 3]));

%% ---- save -----------------------------------------------------------
if ~exist(outdir,'dir'); mkdir(outdir); end
save(fullfile(outdir,'phase9.mat'), '-struct', 'F', '-v7');
fprintf('gen_phase9: wrote %s (%d fields)\n', fullfile(outdir,'phase9.mat'), numel(fieldnames(F)));
end
