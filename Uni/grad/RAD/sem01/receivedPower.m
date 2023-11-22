function Pr = receivedPower(Pt, freq, G, RCS, R)
%
% Inputs:
        % pt        == input peak power in Watts 
        % freq      == radar operating frequency in Hz
        % G         == antenna gain
        % RCS       == radar cross section in meter squared
        % R         == range to target in m
%    
% Outputs:
        % Pr        == received peak power in Watts  
% 
c = 3e8; % speed of light
lambda = c./freq; % wavelength
Pr = Pt.*G.^2.*lambda.^2.*RCS./((4*pi)^3.*R.^4);

end