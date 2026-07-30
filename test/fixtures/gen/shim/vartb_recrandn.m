function out = vartb_recrandn(varargin)
% Recording wrapper around randn, used only by the fixture generator.
% Draws normally, then appends the draw to a global cell array so the exact
% same sequence can be injected into the Julia port (M2). Behaviourally
% identical to randn from the caller's point of view.
global VARTB_DRAWS
out = randn(varargin{:});
if isempty(VARTB_DRAWS); VARTB_DRAWS = {}; end
VARTB_DRAWS{end+1} = out;
end
