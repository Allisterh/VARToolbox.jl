function gen_phase3b(outdir)
% Closes the gap recorded at the end of the Phase 3 log: the whereToStart>1
% branch of SignRestrictions.
%
% That branch is the one that combines an instrument-identified first column
% of B with sign restrictions on the remaining shocks. It is the ONLY
% consumer of the randn(dy,1) draw site, and the only path that reads
% VAR.sigma_b rather than VAR.sigma. The Phase 3 fixtures exercised only
% square SIGN matrices, so none of it was covered.
%
% The two-stage sequence, per the toolbox's documented usage:
%   1. VARir with ident='iv'  -> populates VAR.Biv and VAR.sigma_b
%   2. SR with a SIGN matrix that OMITS that shock's column
%
% Usage:  gen_phase3b('/abs/path/to/test/fixtures')

if ~exist('outdir','var'); error('gen_phase3b: outdir is required'); end
root = '/Users/hiro/GitHub/VAR-Toolbox/v3dot1';
addpath(genpath(root));
shimdir = fullfile(fileparts(mfilename('fullpath')),'shim');
addpath(shimdir);
rehash;

global VARTB_DRAWS

F = struct();

%% ---- Gertler-Karadi VAR with the instrument -------------------------
cd(fullfile(root,'Replic','GK2015'));
[xlsdata, xlstext] = xlsread('GK2015_Data.xlsx','Sheet1');
data = Num2NaN(xlsdata);
vnames = xlstext(2,2:end);
VARvnames = {'gs1';'logcpi';'logip';'ebp';};
ENDO = nan(size(data,1),numel(VARvnames));
for ii=1:numel(VARvnames)
    ENDO(:,ii) = data(:,strcmp(vnames,VARvnames{ii}));
end
IV = data(:,strcmp(vnames,'ff4_tc'));
F.ENDO = ENDO; F.IV = IV;

[VAR, VARopt] = VARmodel(ENDO,12,1);
VAR.IV = IV;
VARopt.nsteps = 12;

% Stage 1: instrument identification populates Biv and sigma_b.
VARopt.ident = 'iv';
[~, VAR] = VARir(VAR,VARopt);
F.Biv     = VAR.Biv;
F.sigma_b = VAR.sigma_b;

%% ---- Stage 2: SignRestrictions on the remaining three shocks --------
% SIGN is 4x3 -- one column short of square, so whereToStart = 1+4-3 = 2
% and the first column is taken as already identified by the instrument.
SIGN = [ 0  0  0
         1  0  0
         1  0  0
         0  1  0];
F.SIGN = SIGN;
F.whereToStart = 1 + size(SIGN,1) - size(SIGN,2);

VARTB_DRAWS = {};
rng(5001);
VARopt.sr_hor = 1;
B = SignRestrictions(SIGN,VAR,VARopt);
F.sr_B = B;
F.sr_draws = VARTB_DRAWS;
F.sr_found = double(~isempty(B));

% The draw sequence here MUST include randn(dy,1) calls -- that is the
% whole point of this fixture. Record the shapes so the test can assert it.
sz = zeros(numel(VARTB_DRAWS),2);
for k = 1:numel(VARTB_DRAWS); sz(k,:) = size(VARTB_DRAWS{k}); end
F.sr_draw_sizes = sz;

%% ---- The dimension checks that guard this branch --------------------
% SIGN with more columns than rows is rejected outright.
F.err_wide_SIGN = 0;
try
    SignRestrictions(ones(3,4),VAR,VARopt);
catch
    F.err_wide_SIGN = 1;
end

% Biv and SIGN must have coherent sizes: with Biv occupying one column,
% SIGN must be exactly one column short of square. A 4x2 SIGN implies
% whereToStart=3, i.e. two already-identified columns, but only one is
% present -> error.
F.err_incoherent = 0;
try
    SignRestrictions([0 0; 1 0; 1 0; 0 1],VAR,VARopt);
catch
    F.err_incoherent = 1;
end

rmpath(shimdir);

%% ---- save -----------------------------------------------------------
if ~exist(outdir,'dir'); mkdir(outdir); end
out = fullfile(outdir,'phase3b.mat');
save(out,'-struct','F','-v7');
fprintf('wrote %s (%d fields)\n', out, numel(fieldnames(F)));
end
