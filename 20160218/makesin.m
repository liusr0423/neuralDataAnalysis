function [m,t] = makesin(sr,f,len)
% MAKECOS Makes a 'virtually' analog sine
%   [m,t] = makesin(f,len)
%   sr: sampling rate in hertz
%   f:   frequency of the sine in hertz
%   len: length in seconds
%   m:   the sine signal
%   t:   the time vector on which m is defined

% Default value
if nargin < 2; len = 2; end


t   = (0:1/sr:(len)-1/sr)'; % time index
m   = sin(2*pi*f*t); % our signal
