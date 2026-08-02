function gen_phase6(outdir)
% Golden-master fixtures for Phase 6: the four published replications.
%
%   Stock & Watson (2001, JEP)      -- Cholesky, IR + VD + bands
%   Blanchard & Quah (1989, AER)    -- long-run restrictions
%   Uhlig (2005, JME)               -- sign restrictions
%   Gertler & Karadi (2015, AEJ)    -- external instrument
%
% This reproduces the NUMERICAL section of each Replic/*/GO_*.m verbatim
% (same data, lags, deterministics, horizons and identification) and dumps
% the results. The plotting sections are omitted -- they are out of scope
% for this port, and three of the four GO scripts die there anyway on the
% missing export_fig.
%
% All shims are active, so every random draw is recorded and the Julia side
% can reproduce each replication exactly rather than distributionally.
%
% Deviation from the published settings, for fixture size only:
%   Uhlig  ndraws 500 -> 100 (recording every posterior draw of a 73x6
%          coefficient matrix for 500+ draws makes the .mat unwieldy).
% Both languages use the reduced count, so the comparison is unaffected;
% only the width of the reported bands differs from the paper.
%
% Revision note: the first version of this generator also silently left
% VARopt.sr_hor and VARopt.pctg at their defaults (1, 95) where
% GO_Uhlig2005.m sets 6 and 68 -- so the Uhlig fixture replicated a variant
% of the paper, not the paper. Caught by the Phase 8 figure comparison;
% both are now set to the GO script's values.
%
% Usage:  gen_phase6('/abs/path/to/test/fixtures')

if ~exist('outdir','var'); error('gen_phase6: outdir is required'); end
root = '/Users/hiro/GitHub/VAR-Toolbox/v3dot1';
addpath(genpath(root));
shimdir = fullfile(fileparts(mfilename('fullpath')),'shim');
addpath(shimdir);
rehash;

global VARTB_URAND VARTB_DRAWS VARTB_POST

F = struct();

%% =====================================================================
%  Stock & Watson (2001) -- GO_SW2001.m
%  =====================================================================
cd(fullfile(root,'Replic','SW2001'));
X = xlsread('SW2001_Data.xlsx','Sheet1');
F.sw_X = X;
[VAR, VARopt] = VARmodel(X,4,1);
VARopt.nsteps = 24;
VARopt.ident  = 'short';
F.sw_nsteps = VARopt.nsteps;
F.sw_ndraws = VARopt.ndraws;
F.sw_pctg   = VARopt.pctg;

[IR, VAR] = VARir(VAR,VARopt);
F.sw_IR = IR;
F.sw_B  = VAR.B;

VARTB_URAND = {}; rng(1001);
[IRinf,IRsup,IRmed,IRbar] = VARirband(VAR,VARopt);
F.sw_irband_draws = VARTB_URAND;
F.sw_IRinf = IRinf; F.sw_IRsup = IRsup;
F.sw_IRmed = IRmed; F.sw_IRbar = IRbar;

[VD, VAR] = VARvd(VAR,VARopt);
F.sw_VD = VD;

VARTB_URAND = {}; rng(1002);
[VDinf,VDsup,VDmed,VDbar] = VARvdband(VAR,VARopt);
F.sw_vdband_draws = VARTB_URAND;
F.sw_VDinf = VDinf; F.sw_VDsup = VDsup;
F.sw_VDmed = VDmed; F.sw_VDbar = VDbar;

% Table 1.B of the paper, exactly as GO_SW2001.m assembles it.
T = zeros(12,3);
T(1:4,  :) = [VD(1,:,1); VD(4,:,1); VD(8,:,1); VD(12,:,1)];
T(5:8,  :) = [VD(1,:,2); VD(4,:,2); VD(8,:,2); VD(12,:,2)];
T(9:12, :) = [VD(1,:,3); VD(4,:,3); VD(8,:,3); VD(12,:,3)];
F.sw_table1B = T;

%% =====================================================================
%  Blanchard & Quah (1989) -- GO_BQ1989.m
%  =====================================================================
cd(fullfile(root,'Replic','BQ1989'));
X = xlsread('BQ1989_Data.xlsx','Sheet1');
F.bq_X = X;
[VAR, VARopt] = VARmodel(X,8,1);
VARopt.nsteps = 40;
VARopt.ident  = 'long';
F.bq_nsteps = VARopt.nsteps;

[IR, VAR] = VARir(VAR,VARopt);
F.bq_IR = IR;
F.bq_B  = VAR.B;

VARTB_URAND = {}; rng(2001);
[IRinf,IRsup,IRmed,IRbar] = VARirband(VAR,VARopt);
F.bq_irband_draws = VARTB_URAND;
F.bq_IRinf = IRinf; F.bq_IRsup = IRsup;
F.bq_IRmed = IRmed; F.bq_IRbar = IRbar;

% Figure 2 of the paper: cumulated supply and demand responses.
F.bq_fig2_supply_gdp = cumsum(IR(:,1,1));
F.bq_fig2_supply_u   = IR(:,2,1);
F.bq_fig2_demand_gdp = cumsum(-IR(:,1,2));
F.bq_fig2_demand_u   = -IR(:,2,2);

%% =====================================================================
%  Gertler & Karadi (2015) -- GO_GK2015.m
%  =====================================================================
cd(fullfile(root,'Replic','GK2015'));
[xlsdata, xlstext] = xlsread('GK2015_Data.xlsx','Sheet1');
data = Num2NaN(xlsdata);
vnames = xlstext(2,2:end);
VARvnames = {'gs1';'logcpi';'logip';'ebp';};
IVvnames  = {'ff4_tc'};
ENDO = nan(size(data,1),numel(VARvnames));
for ii=1:numel(VARvnames)
    ENDO(:,ii) = data(:,strcmp(vnames,VARvnames{ii}));
end
IV = data(:,strcmp(vnames,IVvnames{1}));
F.gk_ENDO = ENDO; F.gk_IV = IV;

[VAR, VARopt] = VARmodel(ENDO,12,1);
VAR.IV = IV;
VARopt.nsteps = 48;
VARopt.ndraws = 200;
VARopt.method = 'wild';
F.gk_nsteps = VARopt.nsteps; F.gk_ndraws = VARopt.ndraws;

% Cholesky block
VARopt.ident = 'short';
[IR, VAR] = VARir(VAR,VARopt);
F.gk_IR_chol = IR;
F.gk_B_chol  = VAR.B;
VARTB_URAND = {}; rng(3001);
[INF,SUP,MED,BAR] = VARirband(VAR,VARopt);
F.gk_chol_draws = VARTB_URAND;
F.gk_chol_INF = INF; F.gk_chol_SUP = SUP;
F.gk_chol_MED = MED; F.gk_chol_BAR = BAR;

% External-instrument block
VARopt.ident = 'iv';
[IRiv, VAR] = VARir(VAR,VARopt);
F.gk_IR_iv = IRiv;
F.gk_B_iv  = VAR.B;
F.gk_Biv   = VAR.Biv;
VARTB_URAND = {}; rng(3002);
[INFiv,SUPiv,MEDiv,BARiv] = VARirband(VAR,VARopt);
F.gk_iv_draws = VARTB_URAND;
F.gk_iv_INF = INFiv; F.gk_iv_SUP = SUPiv;
F.gk_iv_MED = MEDiv; F.gk_iv_BAR = BARiv;

%% =====================================================================
%  Uhlig (2005) -- GO_Uhlig2005.m
%  =====================================================================
cd(fullfile(root,'Replic','Uhlig2005'));
[xlsdata, xlstext] = xlsread('Uhlig2005_Data.xlsx','Sheet1');
data = Num2NaN(xlsdata);
vnames = xlstext(2,2:end);
% GO_Uhlig2005.m rescales five of the six series by 100.
tempnames = {'y','pi','comm','nbres','res'};
for ii=1:numel(tempnames)
    j = strcmp(vnames,tempnames{ii});
    data(:,j) = 100*data(:,j);
end
X = data;
F.uh_X = X;
F.uh_vnames = vnames;

[VAR, VARopt] = VARmodel(X,12,1);
VARopt.nsteps = 60;
VARopt.ndraws = 100;      % reduced from 500 for fixture size; see header
VARopt.sr_hor = 6;        % GO_Uhlig2005.m: restrictions imposed for 6 months
VARopt.pctg   = 68;       % GO_Uhlig2005.m: 68% credible intervals
F.uh_nsteps = VARopt.nsteps; F.uh_ndraws = VARopt.ndraws;
F.uh_sr_mod = VARopt.sr_mod;
F.uh_sr_hor = VARopt.sr_hor; F.uh_pctg = VARopt.pctg;

SIGN = [ 0,0,0,0,0,0
        -1,0,0,0,0,0
        -1,0,0,0,0,0
         0,0,0,0,0,0
        -1,0,0,0,0,0
         1,0,0,0,0,0];
F.uh_SIGN = SIGN;

VARTB_DRAWS = {}; VARTB_POST = {};
rng(4001);
SRout = SR(VAR,SIGN,VARopt);
F.uh_sr_draws  = VARTB_DRAWS;
F.uh_post      = VARTB_POST;
F.uh_Ball      = SRout.Ball;
F.uh_B         = SRout.B;
F.uh_Bmed      = SRout.Bmed;
F.uh_IRmed     = SRout.IRmed;
F.uh_IRinf     = SRout.IRinf;
F.uh_IRsup     = SRout.IRsup;
F.uh_VDmed     = SRout.VDmed;

rmpath(shimdir);

%% ---- save -----------------------------------------------------------
if ~exist(outdir,'dir'); mkdir(outdir); end
out = fullfile(outdir,'phase6.mat');
save(out,'-struct','F','-v7');
d = dir(out);
fprintf('wrote %s (%d fields, %.1f MB)\n', out, numel(fieldnames(F)), d.bytes/1e6);
end
