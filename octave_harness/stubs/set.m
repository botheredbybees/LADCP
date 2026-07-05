function varargout = set(h, varargin)
% octave_harness stub (CONTINUATION_PLAN.md M2): plotraw.m calls
% set(gca,'YTick',...) directly on our fake numeric handle (from the gca.m
% stub), which the real set() rejects. Per the plan's guidance, only
% short-circuit the graphics-handle case; anything else (e.g. set() used on
% a real object/struct) still goes to the real builtin via builtin('set',...).
  if nargin > 0 && isnumeric(h)
    if nargout > 0
      varargout{1} = [];
    end
    return
  end
  [varargout{1:nargout}] = builtin('set', h, varargin{:});
end
