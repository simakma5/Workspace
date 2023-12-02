function out=green( z, x, a, k, d)

% evaluates Green's function for numeric integration

% z - integration variable
% x - distance between the middle of the segmant and the point in which potential computed
% a - radius of the antnna
% k - wavenumber

%R   = sqrt(a^2 + d^2 + (x-z).^2);

if d==0;

    eta  = x-z;
%     G0   = (2./(pi*sqrt(eta.^2)+4*a^2)).*ellipke((4*a^2)./(eta.^2+4*a^2));
    G1   = (exp(-j*k*sqrt(eta.^2+a^2))-1)./sqrt(eta.^2+a^2);
    G0   = 1./sqrt(eta.^2+a^2);
    out   = (G0+G1)./(4*pi);
    %out =  exp(-j*k*sqrt(eta.^2+a^2))./(4*pi*sqrt(eta.^2+a^2));
   
    
else    
 R   = sqrt(a^2 + d^2 + (x-z).^2);
 out = exp(-j*k*R)./(4*pi*R);

end;

