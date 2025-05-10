function y = isoctave()
% Usage: y = isoctave()

y = exist('OCTAVE_VERSION', 'builtin') ~= 0;
