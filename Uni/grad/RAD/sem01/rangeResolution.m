function delta_R = rangeResolution(tau)
% This function computes radar range resolution in meters
%
% Inputs:
    % var can be either
        % var == Bandwidth in Hz
        % var == Pulsewidth in seconds
%    
% Outputs:
    % delta_R == range resolution in meters
%    
% Bandwidth may be equal to (1/pulse width)==> indicator = seconds
%
c = 3e8; % speed of light
delta_R = c*tau/2;

end