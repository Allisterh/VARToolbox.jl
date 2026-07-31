function gen_phase7(outdir)
% Golden-master fixtures for the Stats/ folder.
%
% Deterministic throughout. Input is built from an explicit LCG rather than
% rand, for the same byte-reproducibility reason as gen_phase2.m.
%
% Usage:  gen_phase7('/abs/path/to/test/fixtures')

if ~exist('outdir','var'); error('gen_phase7: outdir is required'); end
root = '/Users/hiro/GitHub/VAR-Toolbox/v3dot1';
addpath(genpath(root));

F = struct();

%% ---- input ----------------------------------------------------------
T = 80;
lcg = zeros(T*4,1); s = 24680;
for i = 1:T*4
    s = mod(1103515245*s + 12345, 2^31);
    lcg(i) = s/2^31 - 0.5;
end
X = [ 10 + reshape(lcg(1:T),T,1)*3, ...
       8 + reshape(lcg(T+1:2*T),T,1)*2, ...
      12 + reshape(lcg(2*T+1:3*T),T,1)*4, ...
       6 + reshape(lcg(3*T+1:4*T),T,1) ];
F.X = X;

% A copy with NaNs, to exercise the NaN-handling differences.
XN = X; XN(5,2) = NaN; XN(20,3) = NaN; XN(40,1) = NaN;
F.XN = XN;

%% ---- moving-window transforms ---------------------------------------
F.movavg_5    = MovAvg(X,5);
F.movavg_12   = MovAvg(X,12);
F.movsum_5    = MovSum(X,5);
F.movstd_5    = MovStd(X,5);
F.movavgcent_4 = MovAvgCent(X,4);
% MovAvgCent uses nanmean, the trailing ones use plain mean -- pin both.
F.movavg_5_nan     = MovAvg(XN,5);
F.movavgcent_4_nan = MovAvgCent(XN,4);

F.movcorr_10_1     = MovCorr(X,10,1);
F.movcorr_10_3     = MovCorr(X,10,3);
F.movcorrcent_8_2  = MovCorrCent(X,8,2);

%% ---- censor ---------------------------------------------------------
F.censor_10 = censor(X,10);
F.censor_25 = censor(X,25);

%% ---- correlation tables ---------------------------------------------
F.corrtable      = CorrTable(X,{'a','b','c','d'});
F.corrtable_pair = CorrTable(XN,{'a','b','c','d'},1);

%% ---- PairCorr -------------------------------------------------------
% Needs strictly positive levels (it takes logs unconditionally).
F.paircorr      = PairCorr(X);
F.paircorr_abs  = PairCorr(X,{'a','b','c','d'},1);

%% ---- CrossCorrelation -----------------------------------------------
[cc, lags] = CrossCorrelation(X(:,1),X(:,2:end),6);
F.crosscorr = cc;
F.crosscorr_lags = lags;

%% ---- error contracts ------------------------------------------------
F.err_window_zero = 0;
try, MovAvg(X,0); catch, F.err_window_zero = 1; end
F.err_window_big = 0;
try, MovAvg(X,1000); catch, F.err_window_big = 1; end

%% ---- save -----------------------------------------------------------
if ~exist(outdir,'dir'); mkdir(outdir); end
out = fullfile(outdir,'phase7.mat');
save(out,'-struct','F','-v7');
fprintf('wrote %s (%d fields)\n', out, numel(fieldnames(F)));
end
