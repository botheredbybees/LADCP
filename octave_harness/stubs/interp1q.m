function yi = interp1q(x, y, xi)
% interp1q -- MATLAB builtin, not implemented in Octave (confirmed under
% gnuoctave/octave:9.2.0). MATLAB's interp1q assumes x is a monotonically
% increasing column vector and does plain linear interpolation with no
% extrapolation (NaN outside range) -- equivalent to Octave's interp1
% with the default 'linear' method.
  yi = interp1(x, y, xi, 'linear');
end
