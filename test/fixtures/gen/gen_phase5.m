function gen_phase5(outdir)
% Golden-master fixtures for Phase 5: LPoption, LPmakex, LPmodel.
%
% Fully deterministic -- local projections involve no RNG.
%
% Usage:  gen_phase5('/abs/path/to/test/fixtures')

if ~exist('outdir','var'); error('gen_phase5: outdir is required'); end
root = '/Users/hiro/GitHub/VAR-Toolbox/v3dot1';
addpath(genpath(root));

F = struct();

%% ---- LPoption defaults ----------------------------------------------
o = LPoption;
F.opt_nsteps = o.nsteps;
F.opt_pctg   = o.pctg;
F.opt_impact = o.impact;
F.opt_EXOG_empty = double(isempty(o.EXOG));
F.opt_LAGS_empty = double(isempty(o.LAGS));
% LPoption.m defines NO nlag / const fields, so LPmodel errors on an
% unmodified options struct. Pinned as a contract.
F.opt_has_nlag  = double(isfield(o,'nlag'));
F.opt_has_const = double(isfield(o,'const'));

%% ---- data -----------------------------------------------------------
cd(fullfile(root,'Replic','SW2001'));
SW = xlsread('SW2001_Data.xlsx','Sheet1');
F.SW = SW;

%% ---- LPmakex --------------------------------------------------------
A = SW(1:20,:);
F.A = A;
F.lpmakex_1 = LPmakex(A,1);
F.lpmakex_3 = LPmakex(A,3);
F.lpmakex_1col = LPmakex(A(:,1),2);

%% ---- LPmodel: with LAGS --------------------------------------------
ENDO = SW(:,1);       % inflation
S    = SW(:,3);       % fed funds as the "shock" series
LAGS = SW;
opt = LPoption;
opt.nsteps = 16;
opt.nlag   = 4;
opt.const  = 1;
opt.LAGS   = LAGS;
opt.impact = 0;
[IR,INF,SUP,LPout] = LPmodel(ENDO,S,opt);
F.lp1_IR = IR; F.lp1_INF = INF; F.lp1_SUP = SUP;
F.lp1_nvar = LPout.nvar; F.lp1_ntotcoeff = LPout.ntotcoeff;
F.lp1_nsteps = opt.nsteps; F.lp1_nlag = opt.nlag;
for hh = 1:opt.nsteps
    a = LPout.(sprintf('h%d',hh));
    q = sprintf('lp1_h%d_', hh);
    F.([q 'nobs'])    = a.nobs;
    F.([q 'beta'])    = a.beta;
    F.([q 'tstat'])   = a.tstat;
    F.([q 'bstd'])    = a.bstd;
    F.([q 'bstd_HW']) = a.bstd_HW;   % NB: holds NW values (overwritten upstream)
    F.([q 'tprob'])   = a.tprob;
    F.([q 'rsqr'])    = a.rsqr;
    F.([q 'dw'])      = a.dw;
    F.([q 'sigma'])   = a.sigma;
end

%% ---- LPmodel: unit shock -------------------------------------------
opt2 = opt; opt2.impact = 1;
[IR2,INF2,SUP2,~] = LPmodel(ENDO,S,opt2);
F.lp2_IR = IR2; F.lp2_INF = INF2; F.lp2_SUP = SUP2;

%% ---- LPmodel: no LAGS (=> no intercept at all) ---------------------
opt3 = LPoption;
opt3.nsteps = 8; opt3.nlag = 4; opt3.const = 1;
opt3.LAGS = []; opt3.impact = 0;
[IR3,INF3,SUP3,LPout3] = LPmodel(ENDO,S,opt3);
F.lp3_IR = IR3; F.lp3_INF = INF3; F.lp3_SUP = SUP3;
F.lp3_ntotcoeff = LPout3.ntotcoeff;
F.lp3_h1_beta = LPout3.h1.beta;
F.lp3_h1_nobs = LPout3.h1.nobs;
F.lp3_nsteps = opt3.nsteps;

%% ---- LPmodel: with EXOG --------------------------------------------
T = size(SW,1);
lcg = zeros(T,1); s0 = 777;
for i = 1:T; s0 = mod(1103515245*s0 + 12345, 2^31); lcg(i) = s0/2^31 - 0.5; end
F.EX = lcg;
opt4 = opt; opt4.EXOG = lcg; opt4.nsteps = 8;
[IR4,INF4,SUP4,LPout4] = LPmodel(ENDO,S,opt4);
F.lp4_IR = IR4; F.lp4_INF = INF4; F.lp4_SUP = SUP4;
F.lp4_ntotcoeff = LPout4.ntotcoeff;
F.lp4_h1_beta = LPout4.h1.beta;
F.lp4_nsteps = opt4.nsteps;

%% ---- LPmodel: different coverage -----------------------------------
opt5 = opt; opt5.pctg = 68; opt5.nsteps = 8;
[IR5,INF5,SUP5,~] = LPmodel(ENDO,S,opt5);
F.lp5_IR = IR5; F.lp5_INF = INF5; F.lp5_SUP = SUP5;
F.lp5_nsteps = opt5.nsteps;

%% ---- error contracts ------------------------------------------------
F.lp_bad_S_errors = 0;
try, LPmodel(ENDO,S(1:end-5),opt); catch, F.lp_bad_S_errors = 1; end
F.lp_bad_impact_errors = 0;
opt6 = opt; opt6.impact = 7;
try, LPmodel(ENDO,S,opt6); catch, F.lp_bad_impact_errors = 1; end
F.lp_no_nlag_errors = 0;
opt7 = LPoption; opt7.LAGS = LAGS;
try, LPmodel(ENDO,S,opt7); catch, F.lp_no_nlag_errors = 1; end

%% ---- save -----------------------------------------------------------
if ~exist(outdir,'dir'); mkdir(outdir); end
out = fullfile(outdir,'phase5.mat');
save(out,'-struct','F','-v7');
fprintf('wrote %s (%d fields)\n', out, numel(fieldnames(F)));
end
