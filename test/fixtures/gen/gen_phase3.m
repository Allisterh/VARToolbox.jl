function gen_phase3(outdir)
% Golden-master fixtures for Phase 3: VARir (all five ident schemes),
% VARvd, VARhd, SignRestrictions, SR, VARdrawpost.
%
% M2 in practice: the two randn() call sites inside SignRestrictions are
% routed through a recording wrapper (see gen/shim/), so the exact draw
% sequence MATLAB used is saved alongside the results. Replaying that
% sequence in Julia turns the whole sign-restriction pipeline into an exact
% golden-master check instead of a distributional comparison.
%
% Usage:  gen_phase3('/abs/path/to/test/fixtures')

if ~exist('outdir','var'); error('gen_phase3: outdir is required'); end
root = '/Users/hiro/GitHub/VAR-Toolbox/v3dot1';
addpath(genpath(root));

global VARTB_DRAWS

F = struct();

%% ---- data -----------------------------------------------------------
cd(fullfile(root,'Replic','SW2001'));
SW = xlsread('SW2001_Data.xlsx','Sheet1');
F.SW = SW;

[V, o] = VARmodel(SW,4,1);
o.nsteps = 12;

%% ---- VARir: short / long --------------------------------------------
for nm = {'short','long'}
    oo = o; oo.ident = nm{1};
    [IR, Vout] = VARir(V,oo);
    F.(['ir_' nm{1}])       = IR;
    F.(['B_'  nm{1}])       = Vout.B;
    F.(['PSI_' nm{1}])      = Vout.PSI;
    F.(['Fp_'  nm{1}])      = Vout.Fp;
end

% recurs = 'comp' should agree with 'wold'; both saved
oo = o; oo.ident = 'short'; oo.recurs = 'comp';
F.ir_short_comp = VARir(V,oo);

% impact = 1 (unit shock)
oo = o; oo.ident = 'short'; oo.impact = 1;
F.ir_short_unit = VARir(V,oo);

% shut
oo = o; oo.ident = 'short'; oo.shut = 2;
F.ir_short_shut2 = VARir(V,oo);

%% ---- VARvd / VARhd: short / long / sign -----------------------------
for nm = {'short','long'}
    oo = o; oo.ident = nm{1};
    F.(['vd_' nm{1}]) = VARvd(V,oo);
    HD = VARhd(V,oo);
    F.(['hd_' nm{1} '_shock'])  = HD.shock;
    F.(['hd_' nm{1} '_init'])   = HD.init;
    F.(['hd_' nm{1} '_const'])  = HD.const;
    F.(['hd_' nm{1} '_trend'])  = HD.trend;
    F.(['hd_' nm{1} '_trend2']) = HD.trend2;
    F.(['hd_' nm{1} '_endo'])   = HD.endo;
end

% VD sums to 100 across its SECOND index (shocks), confirming that VD is
% [nsteps x nshocks x nvar] while IR is [nsteps x nvar x nshocks].
F.vd_sum_dim2 = sum(F.vd_short(1,:,1));
F.vd_sum_dim3 = sum(F.vd_short(1,1,:));

% Which ident values each function rejects (contract, not a number).
F.vd_rejects_iv = 0;  F.vd_rejects_exog = 0;
F.hd_rejects_iv = 0;  F.hd_rejects_exog = 0;
oo = o; oo.ident = 'iv';   try, VARvd(V,oo); catch, F.vd_rejects_iv = 1;   end
oo = o; oo.ident = 'exog'; try, VARvd(V,oo); catch, F.vd_rejects_exog = 1; end
oo = o; oo.ident = 'iv';   try, VARhd(V,oo); catch, F.hd_rejects_iv = 1;   end
oo = o; oo.ident = 'exog'; try, VARhd(V,oo); catch, F.hd_rejects_exog = 1; end

% VARir with ident='sign' but no B must raise.
F.ir_sign_no_B_errors = 0;
oo = o; oo.ident = 'sign';
try, VARir(V,oo); catch, F.ir_sign_no_B_errors = 1; end

%% ---- VARir with ident='exog' ----------------------------------------
% Needs an exogenous variable. Use a deterministic LCG series.
T = size(SW,1);
lcg = zeros(T,1); s = 999;
for i = 1:T; s = mod(1103515245*s + 12345, 2^31); lcg(i) = s/2^31 - 0.5; end
F.EX = lcg;
[Vx, ox] = VARmodel(SW,4,1,lcg,0);
ox.nsteps = 12; ox.ident = 'exog';
[IRx, Vxo] = VARir(Vx,ox);
F.ir_exog = IRx;
F.B_exog  = Vxo.B;
F.vx_ntotcoeff = Vx.ntotcoeff;
F.vx_nvar_ex   = Vx.nvar_ex;

% HD with an exogenous variable (exercises the HDexo branch)
oh = ox; oh.ident = 'short';
HDx = VARhd(Vx,oh);
F.hdx_exo   = HDx.exo;
F.hdx_endo  = HDx.endo;
F.hdx_shock = HDx.shock;

%% ---- VARir with ident='iv' (Gertler-Karadi) -------------------------
cd(fullfile(root,'Replic','GK2015'));
[gkdata, gktext] = xlsread('GK2015_Data.xlsx','Sheet1');
% GO_GK2015.m applies Num2NaN BEFORE selecting columns -- the sheet encodes
% missing values as the sentinel 123456789, and the ff4_tc instrument is
% only populated over part of the sample. Skipping this step leaves sentinel
% values in the instrument and makes the IV first stage numerically
% singular (RCOND ~1e-16), which is how this was caught.
gkdata = Num2NaN(gkdata);
F.gk_raw = gkdata;
vnames = gktext(2,2:end);
idxIV  = find(strcmp(vnames,'ff4_tc'));
idxEND = [find(strcmp(vnames,'gs1')), find(strcmp(vnames,'logcpi')), ...
          find(strcmp(vnames,'logip')), find(strcmp(vnames,'ebp'))];
F.gk_idxIV = idxIV; F.gk_idxEND = idxEND;
ENDO = gkdata(:,idxEND);
IVv  = gkdata(:,idxIV);
if any(isnan(ENDO(:)))
    error('gen_phase3: ENDO contains NaN after Num2NaN -- check column selection');
end
[Vg, og] = VARmodel(ENDO,12,1);
Vg.IV = IVv;
og.nsteps = 12; og.ident = 'iv';
[IRiv, Vgo] = VARir(Vg,og);
F.gk_ENDO = ENDO;  F.gk_IV = IVv;
F.ir_iv    = IRiv;
F.B_iv     = Vgo.B;
F.Biv      = Vgo.Biv;
F.sigma_b  = Vgo.sigma_b;
F.gk_fs_beta  = Vgo.FirstStage.beta;
F.gk_fs_yhat  = Vgo.FirstStage.yhat;
F.gk_fs_rsqr  = Vgo.FirstStage.rsqr;
F.gk_sigma = Vg.sigma;
F.gk_Ft    = Vg.Ft;

% Order matters: swapping the instrumented variable out of position 1 must
% change the answer. Pinned so the convention is provably load-bearing.
ENDO2 = ENDO(:,[2 1 3 4]);
[Vg2, og2] = VARmodel(ENDO2,12,1);
Vg2.IV = IVv; og2.nsteps = 12; og2.ident = 'iv';
F.ir_iv_swapped = VARir(Vg2,og2);

%% ---- SignRestrictions with recorded draws (M2) ----------------------
% Shadow the upstream SignRestrictions with the recording shim.
shimdir = fullfile(fileparts(mfilename('fullpath')),'shim');
addpath(shimdir);                     % prepends -> shadows v3dot1 copy
rehash;

cd(fullfile(root,'Replic','SW2001'));
SIGN = [ 1  0  0
        -1  0  0
         1  0  0];
F.SIGN = SIGN;

VARTB_DRAWS = {};
rng(42);
osr = o; osr.ident = 'sign'; osr.sr_hor = 1;
Bsr = SignRestrictions(SIGN,V,osr);
F.sr_B = Bsr;
F.sr_draws = VARTB_DRAWS;
F.sr_ndraws_used = numel(VARTB_DRAWS);

% sr_hor > 1 exercises the VARir-inside-the-loop path
VARTB_DRAWS = {};
rng(43);
osr2 = o; osr2.ident = 'sign'; osr2.sr_hor = 4;
Bsr2 = SignRestrictions(SIGN,V,osr2);
F.sr2_B = Bsr2;
F.sr2_draws = VARTB_DRAWS;

%% ---- SR end to end, sr_mod = 0 (model fixed) ------------------------
% sr_mod=0 avoids VARdrawpost, so the ONLY randomness is the recorded
% rotation draws -> the entire SR pipeline is exactly reproducible.
VARTB_DRAWS = {};
rng(7);
osr3 = o; osr3.ndraws = 20; osr3.sr_mod = 0; osr3.nsteps = 8;
SRout = SR(V,SIGN,osr3);
F.SR_draws  = VARTB_DRAWS;
F.SR_ndraws = osr3.ndraws;
F.SR_nsteps = osr3.nsteps;
F.SR_Ball   = SRout.Ball;
F.SR_Bmed   = SRout.Bmed;
F.SR_B      = SRout.B;
F.SR_IRall  = SRout.IRall;
F.SR_VDall  = SRout.VDall;
F.SR_IRmed  = SRout.IRmed;
F.SR_IRinf  = SRout.IRinf;
F.SR_IRsup  = SRout.IRsup;
F.SR_VDmed  = SRout.VDmed;
F.SR_VDinf  = SRout.VDinf;
F.SR_VDsup  = SRout.VDsup;
F.SR_IR     = SRout.IR;
F.SR_VD     = SRout.VD;
F.SR_HD_shock = SRout.HD.shock;
F.SR_HD_endo  = SRout.HD.endo;
F.SR_HDall_endo = SRout.HDall.endo;

% Pin the off-by-one: which accepted draw does SRout.B equal, and which
% does SRout.IR come from? (See PORTING_LOG / SR docstring.)
selB = NaN;
for k = 1:osr3.ndraws
    if isequal(SRout.Ball(:,:,k), SRout.B); selB = k; break; end
end
osel = osr3; osel.ident = 'sign';
bestIR = NaN; bestd = Inf;
for k = 1:osr3.ndraws
    Vk = V; Vk.B = SRout.Ball(:,:,k);
    IRk = VARir(Vk,osel);
    d = max(abs(SRout.IR(:)-IRk(:)));
    if d < bestd; bestd = d; bestIR = k; end
end
F.SR_selB      = selB;
F.SR_selIR     = bestIR;
F.SR_selIRresid = bestd;

rmpath(shimdir);

%% ---- VARdrawpost ----------------------------------------------------
% The draw itself uses wishrnd/mvnrnd (Statistics Toolbox). Save one draw
% so the DETERMINISTIC part downstream of it -- the reshape, F_draw and
% Fcomp_draw -- can be validated exactly by injection.
rng(11);
[sigma_draw, Ft_draw, F_draw, Fcomp_draw] = VARdrawpost(V);
F.vdp_sigma_draw = sigma_draw;
F.vdp_Ft_draw    = Ft_draw;
F.vdp_F_draw     = F_draw;
F.vdp_Fcomp_draw = Fcomp_draw;

%% ---- save -----------------------------------------------------------
if ~exist(outdir,'dir'); mkdir(outdir); end
out = fullfile(outdir,'phase3.mat');
save(out,'-struct','F','-v7');
fprintf('wrote %s (%d fields)\n', out, numel(fieldnames(F)));
end
