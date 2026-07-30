function gen_phase4(outdir)
% Golden-master fixtures for Phase 4: VARirband, VARvdband.
%
% Uses the recording shims (gen/shim/) so the exact rand() sequence MATLAB
% consumed is saved. Replaying it in Julia makes the bootstrap bands an
% exact golden-master check rather than a distributional comparison -- the
% whole point of M2.
%
% Usage:  gen_phase4('/abs/path/to/test/fixtures')

if ~exist('outdir','var'); error('gen_phase4: outdir is required'); end
root = '/Users/hiro/GitHub/VAR-Toolbox/v3dot1';
addpath(genpath(root));
shimdir = fullfile(fileparts(mfilename('fullpath')),'shim');
addpath(shimdir);   % prepended -> shadows the v3dot1 copies
rehash;

global VARTB_URAND

F = struct();

%% ---- data -----------------------------------------------------------
cd(fullfile(root,'Replic','SW2001'));
SW = xlsread('SW2001_Data.xlsx','Sheet1');
F.SW = SW;
[V, o] = VARmodel(SW,4,1);
o.nsteps = 12;
o.ndraws = 50;
o.pctg   = 95;
F.nsteps = o.nsteps;  F.ndraws = o.ndraws;  F.pctg = o.pctg;

%% ---- VARirband: standard bootstrap ----------------------------------
VARTB_URAND = {};
rng(101);
oo = o; oo.ident = 'short'; oo.method = 'bs';
[INF,SUP,MED,BAR] = VARirband(V,oo);
F.irb_bs_draws = VARTB_URAND;
F.irb_bs_INF = INF; F.irb_bs_SUP = SUP;
F.irb_bs_MED = MED; F.irb_bs_BAR = BAR;

%% ---- VARirband: wild bootstrap --------------------------------------
VARTB_URAND = {};
rng(102);
oo = o; oo.ident = 'short'; oo.method = 'wild';
[INF,SUP,MED,BAR] = VARirband(V,oo);
F.irb_wild_draws = VARTB_URAND;
F.irb_wild_INF = INF; F.irb_wild_SUP = SUP;
F.irb_wild_MED = MED; F.irb_wild_BAR = BAR;

%% ---- VARirband: long-run identification -----------------------------
VARTB_URAND = {};
rng(103);
oo = o; oo.ident = 'long'; oo.method = 'bs';
[INF,SUP,MED,BAR] = VARirband(V,oo);
F.irb_long_draws = VARTB_URAND;
F.irb_long_INF = INF; F.irb_long_SUP = SUP;
F.irb_long_MED = MED; F.irb_long_BAR = BAR;

%% ---- VARirband: different coverage ----------------------------------
VARTB_URAND = {};
rng(104);
oo = o; oo.ident = 'short'; oo.method = 'bs'; oo.pctg = 68;
[INF,SUP,MED,BAR] = VARirband(V,oo);
F.irb_p68_draws = VARTB_URAND;
F.irb_p68_INF = INF; F.irb_p68_SUP = SUP;

%% ---- VARvdband ------------------------------------------------------
VARTB_URAND = {};
rng(201);
oo = o; oo.ident = 'short'; oo.method = 'bs';
[INF,SUP,MED,BAR] = VARvdband(V,oo);
F.vdb_bs_draws = VARTB_URAND;
F.vdb_bs_INF = INF; F.vdb_bs_SUP = SUP;
F.vdb_bs_MED = MED; F.vdb_bs_BAR = BAR;

VARTB_URAND = {};
rng(202);
oo = o; oo.ident = 'short'; oo.method = 'wild';
[INF,SUP,MED,BAR] = VARvdband(V,oo);
F.vdb_wild_draws = VARTB_URAND;
F.vdb_wild_INF = INF; F.vdb_wild_SUP = SUP;
F.vdb_wild_MED = MED; F.vdb_wild_BAR = BAR;

%% ---- VARirband with EXOG -------------------------------------------
T = size(SW,1);
lcg = zeros(T,1); s = 555;
for i = 1:T; s = mod(1103515245*s + 12345, 2^31); lcg(i) = s/2^31 - 0.5; end
F.EX = lcg;
[Vx, ox] = VARmodel(SW,4,1,lcg,0);
ox.nsteps = 8; ox.ndraws = 20; ox.ident = 'short'; ox.method = 'bs';
F.exog_nsteps = ox.nsteps; F.exog_ndraws = ox.ndraws;
VARTB_URAND = {};
rng(301);
[INF,SUP,MED,BAR] = VARirband(Vx,ox);
F.irb_exog_draws = VARTB_URAND;
F.irb_exog_INF = INF; F.irb_exog_SUP = SUP;
F.irb_exog_MED = MED; F.irb_exog_BAR = BAR;

%% ---- VARirband with iv + wild (the branch VARvdband lacks) ----------
cd(fullfile(root,'Replic','GK2015'));
[gkdata, gktext] = xlsread('GK2015_Data.xlsx','Sheet1');
gkdata = Num2NaN(gkdata);
vnames = gktext(2,2:end);
idxIV  = find(strcmp(vnames,'ff4_tc'));
idxEND = [find(strcmp(vnames,'gs1')), find(strcmp(vnames,'logcpi')), ...
          find(strcmp(vnames,'logip')), find(strcmp(vnames,'ebp'))];
ENDO = gkdata(:,idxEND);
IVv  = gkdata(:,idxIV);
[Vg, og] = VARmodel(ENDO,12,1);
Vg.IV = IVv;
og.nsteps = 8; og.ndraws = 10; og.ident = 'iv'; og.method = 'wild';
F.iv_nsteps = og.nsteps; F.iv_ndraws = og.ndraws;
F.gk_ENDO = ENDO; F.gk_IV = IVv;
VARTB_URAND = {};
rng(401);
[INF,SUP,MED,BAR] = VARirband(Vg,og);
F.irb_iv_draws = VARTB_URAND;
F.irb_iv_INF = INF; F.irb_iv_SUP = SUP;
F.irb_iv_MED = MED; F.irb_iv_BAR = BAR;

%% ---- unknown method must raise --------------------------------------
F.band_bad_method_errors = 0;
oo = o; oo.method = 'nonsense';
try, VARirband(V,oo); catch, F.band_bad_method_errors = 1; end

rmpath(shimdir);

%% ---- save -----------------------------------------------------------
if ~exist(outdir,'dir'); mkdir(outdir); end
out = fullfile(outdir,'phase4.mat');
save(out,'-struct','F','-v7');
fprintf('wrote %s (%d fields)\n', out, numel(fieldnames(F)));
end
