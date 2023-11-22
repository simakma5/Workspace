function [fd, tdr] = dopplerFreq(freq, ang, tv)
% This function computes Doppler frequency and time dilation factor ratio
% tau_prime / tau
%
% Inputs:
% freq  == radar operating frequency in Hz
% ang   == target aspect angle in degrees
% tv    == target velocity in m/sec
%
% Outputs:
% fd    == Doppler frequency in Hz
% tdr   == time dilation factor; unitless
%

c = 3e8;
ang_rad = ang*pi./180.;
lambda = c./freq;

fd = 2*tv.*cos(ang_rad)./lambda;
tdr = (c - tv)./(c + tv);

end