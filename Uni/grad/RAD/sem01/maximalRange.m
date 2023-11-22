function Rmax = maximalRange(Pt, G, freq, RCS, PrMin)
%
% Inputs:
        % Pt        == input peak power in Watts 
        % freq      == radar operating frequency in Hz
        % G         == antenna gain
        % RCS       == radar cross section in meter squared
        % PrMin     == minimal received peak power in Watts  
%    
% Outputs:
         % Rmax     == maximal range to target in m

c = 3e8;
lambda = c./freq;
Rmax = nthroot(Pt.*G.^2.*lambda.^2.*RCS./((4*pi)^2.*PrMin), 4);

end