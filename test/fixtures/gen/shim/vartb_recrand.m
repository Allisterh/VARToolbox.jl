function out = vartb_recrand(varargin)
% Recording wrapper around rand, used only by the fixture generator.
% Draws normally, then appends the draw to a global cell array so the exact
% same sequence can be injected into the Julia port (M2).
global VARTB_URAND
out = rand(varargin{:});
if isempty(VARTB_URAND); VARTB_URAND = {}; end
VARTB_URAND{end+1} = out;
end
