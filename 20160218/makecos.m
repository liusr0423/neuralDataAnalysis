function [m,t] = makecos(f,len)
% MAKECOS Makes a 'virtually' analog cosine
%   [m,t] = makecos(f,len)
%
%   f:   frequency of the cosine in hertz
%   len: length in seconds
%   m:   the cosine signal
%   t:   the time vector on which m is defined

% Default value
if nargin < 2; len = 2; end

sr  = 5000; % ignore this
t   = ((-len/2):1/sr:(len/2)).'; % time index
m   = cos(2*pi*f*t); % our signal
