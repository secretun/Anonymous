% gauss() - return a Gaussian window
%
% Usage:
%          >> outvec = gauss(frames,sds);
%
%   frames = window length
%   steep  = steepness (~0+ -> flat; >>10 -> spike)

function outvec = gauss(frames,sds)

outvec = [];
if nargin < 2
  help gauss
  return
end
if sds <=0 | frames < 1
  help gauss
  return
end

incr = 2*sds/(frames-1);
outvec = exp(-(-sds:incr:sds).^2);
