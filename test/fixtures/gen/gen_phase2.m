function gen_phase2(outdir)
% Golden-master fixtures for Phase 2: VARmodel, VARoption, ARDLmodel.
%
% Exercises both a small synthetic matrix (so failures are debuggable) and
% the real Stock-Watson 2001 data at the exact configuration GO_SW2001.m
% uses, so that a Phase 6 end-to-end failure can be localised.
%
% Usage:  gen_phase2('/abs/path/to/test/fixtures')

if ~exist('outdir','var'); error('gen_phase2: outdir is required'); end
root = '/Users/hiro/GitHub/VAR-Toolbox/v3dot1';
addpath(genpath(root));

F = struct();

%% ---- synthetic input -------------------------------------------------
% Deterministic, 60 observations, built from an explicit linear congruential
% generator rather than rand/randn. Two reasons:
%
%   * reproducible without depending on MATLAB's RNG state or algorithm, so
%     regenerating this fixture on another machine gives identical bytes;
%   * white-noise-like, which keeps the regressor matrix well conditioned.
%     A smooth deterministic series (trig, polynomial) makes a variable
%     nearly a linear combination of its own lags -- the first draft of this
%     generator used sin/cos and produced RCOND ~1e-19 designs, against
%     which an rtol=1e-10 comparison is meaningless.
%
% Long enough that nlag=3 with const=3 still leaves positive degrees of
% freedom: VARmodel computes t-probabilities on nobse-ncoeff dof and
% tdis_prb rejects dof<=0.
T = 60;
lcg = zeros(T*5,1);
s = 12345;
for i = 1:T*5
    s = mod(1103515245*s + 12345, 2^31);
    lcg(i) = s / 2^31 - 0.5;          % ~U(-0.5, 0.5)
end
A = [ 2.0 + reshape(lcg(1:T),T,1), ...
     -1.0 + reshape(lcg(T+1:2*T),T,1), ...
      0.5 + reshape(lcg(2*T+1:3*T),T,1) ];
F.A = A;
E = [ reshape(lcg(3*T+1:4*T),T,1), reshape(lcg(4*T+1:5*T),T,1) ];
F.E = E;

%% ---- VARoption defaults ---------------------------------------------
o = VARoption;
F.opt_nsteps = o.nsteps;   F.opt_ndraws = o.ndraws;  F.opt_pctg    = o.pctg;
F.opt_sr_hor = o.sr_hor;   F.opt_sr_rot = o.sr_rot;  F.opt_sr_draw = o.sr_draw;
F.opt_sr_mod = o.sr_mod;   F.opt_mult   = o.mult;    F.opt_quality = o.quality;
F.opt_impact = o.impact;   F.opt_shut   = o.shut;    F.opt_pick    = o.pick;
F.opt_datestype = o.datestype;  F.opt_suptitle = o.suptitle;
F.opt_ident = o.ident;     F.opt_recurs = o.recurs;  F.opt_method  = o.method;
F.opt_frequency = o.frequency;  F.opt_FigSize = o.FigSize;

%% ---- VARmodel: synthetic, no EXOG ------------------------------------
cfg = {2,1; 1,1; 3,1; 2,0; 2,2; 2,3};
for k = 1:size(cfg,1)
    nlag = cfg{k,1}; cst = cfg{k,2};
    V = VARmodel(A,nlag,cst);
    p = sprintf('v%d_', k);
    F.([p 'nlag']) = V.nlag;   F.([p 'const'])     = V.const;
    F.([p 'nobs']) = V.nobs;   F.([p 'nvar'])      = V.nvar;
    F.([p 'ncoeff']) = V.ncoeff; F.([p 'ntotcoeff']) = V.ntotcoeff;
    F.([p 'Ft'])    = V.Ft;    F.([p 'F'])      = V.F;
    F.([p 'sigma']) = V.sigma; F.([p 'resid'])  = V.resid;
    F.([p 'X'])     = V.X;     F.([p 'Y'])      = V.Y;
    F.([p 'Fcomp']) = V.Fcomp; F.([p 'maxEig']) = V.maxEig;
    for j = 1:V.nvar
        e = eval(sprintf('V.eq%d', j));
        q = sprintf('%seq%d_', p, j);
        F.([q 'beta'])  = e.beta;   F.([q 'tstat']) = e.tstat;
        F.([q 'bstd'])  = e.bstd;   F.([q 'tprob']) = e.tprob;
        F.([q 'resid']) = e.resid;  F.([q 'yhat'])  = e.yhat;
        F.([q 'y'])     = e.y;      F.([q 'rsqr'])  = e.rsqr;
        F.([q 'rbar'])  = e.rbar;   F.([q 'sige'])  = e.sige;
        F.([q 'dw'])    = e.dw;
    end
end
F.nconfig = size(cfg,1);

%% ---- VARmodel with EXOG: all three nlag vs nlag_ex branches ----------
ex = {2,2; 3,1; 1,3};    % equal, nlag>nlag_ex, nlag<nlag_ex
for k = 1:size(ex,1)
    nlag = ex{k,1}; nlx = ex{k,2};
    V = VARmodel(A,nlag,1,E,nlx);
    p = sprintf('x%d_', k);
    F.([p 'nlag'])  = V.nlag;    F.([p 'nlag_ex']) = V.nlag_ex;
    F.([p 'nobs'])  = V.nobs;    F.([p 'nvar_ex']) = V.nvar_ex;
    F.([p 'ncoeff']) = V.ncoeff; F.([p 'ntotcoeff']) = V.ntotcoeff;
    F.([p 'Ft'])    = V.Ft;      F.([p 'F'])     = V.F;
    F.([p 'sigma']) = V.sigma;   F.([p 'resid']) = V.resid;
    F.([p 'X'])     = V.X;       F.([p 'Y'])     = V.Y;
    F.([p 'X_EX'])  = V.X_EX;    F.([p 'Fcomp']) = V.Fcomp;
    F.([p 'maxEig']) = V.maxEig;
end
F.nexog = size(ex,1);

%% ---- VARmodel on real data (SW2001, exactly as GO_SW2001.m) ---------
cd(fullfile(root,'Replic','SW2001'));
[xlsdata, ~] = xlsread('SW2001_Data.xlsx','Sheet1');
SW = xlsdata;
F.SW = SW;
V = VARmodel(SW,4,1);
F.sw_nobs = V.nobs;  F.sw_nvar = V.nvar;
F.sw_ntotcoeff = V.ntotcoeff;
F.sw_Ft = V.Ft;      F.sw_F = V.F;
F.sw_sigma = V.sigma;  F.sw_resid = V.resid;
F.sw_Fcomp = V.Fcomp;  F.sw_maxEig = V.maxEig;
F.sw_X = V.X;          F.sw_Y = V.Y;
for j = 1:V.nvar
    e = eval(sprintf('V.eq%d', j));
    q = sprintf('sw_eq%d_', j);
    F.([q 'beta']) = e.beta;  F.([q 'tstat']) = e.tstat;
    F.([q 'tprob']) = e.tprob; F.([q 'rsqr']) = e.rsqr;
    F.([q 'dw']) = e.dw;      F.([q 'sige']) = e.sige;
end

%% ---- ARDLmodel ------------------------------------------------------
y1 = A(:,1);
D = ARDLmodel(y1,2,1);
F.ardl_nobs = D.nobs;  F.ardl_nvar = D.nvar;  F.ardl_ncoeff = D.ncoeff;
F.ardl_beta = D.beta;  F.ardl_yhat = D.yhat;  F.ardl_resid = D.resid;
F.ardl_sige = D.sige;  F.ardl_bstd = D.bstd;  F.ardl_bint = D.bint;
F.ardl_tstat = D.tstat; F.ardl_tprob = D.tprob;
F.ardl_rsqr = D.rsqr;  F.ardl_rbar = D.rbar;  F.ardl_dw = D.dw;
F.ardl_F = D.F;        F.ardl_sigbeta = D.sigbeta;
F.ardl_x = D.x;        F.ardl_y = D.y;

D2 = ARDLmodel(y1,2,1,E,1);
F.ardl2_nobs = D2.nobs;  F.ardl2_nvar = D2.nvar;
F.ardl2_beta = D2.beta;  F.ardl2_resid = D2.resid;
F.ardl2_tstat = D2.tstat; F.ardl2_tprob = D2.tprob;
F.ardl2_rsqr = D2.rsqr;  F.ardl2_dw = D2.dw;
F.ardl2_x = D2.x;        F.ardl2_y = D2.y;

D3 = ARDLmodel(y1,2,0);
F.ardl3_beta = D3.beta;  F.ardl3_nvar = D3.nvar;  F.ardl3_rsqr = D3.rsqr;

%% ---- save -----------------------------------------------------------
if ~exist(outdir,'dir'); mkdir(outdir); end
out = fullfile(outdir,'phase2.mat');
save(out,'-struct','F','-v7');
fprintf('wrote %s (%d fields)\n', out, numel(fieldnames(F)));
end
