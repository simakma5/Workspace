function dy=derivace(x,y)
% DERIVACE - numericka derivace funkce dane tabulkou
% dy=derivace(x,y)
% x ... vektor hodnot nezavisle promenne, ekvidistantni (delka N)
% y ... vektor funkcnich hodnot (delka N)
% dy ... vektor pribliznych hodnot derivace v bodech [xi,yi] (delka N)

N = length(x);
if N~=length(y)
    error('Vstupni vektory nejsou stejne dlouhe!')
end
if N<3
    error('Zadali jste malo bodu - potrebuji alespon TRI hodnoty!')
end
h = (x(2)-x(1)); % pro ekvidistantn� krok
%...jinak TEST, zda je krok vsude stejny!...

dy = y;  % priprava, 'dy' je stejne velke jako 'y'
% zmena hodnot 'dy' na spravne:
dy(1) = (-3*y(1)+4*y(2)-y(3))/(2*h); % pocatek
for i=2:N-1
    dy(i) = (y(i+1)-y(i-1))/(2*h); % vnitrni body
end
dy(N) = (y(N-2)-4*y(N-1)+3*y(N))/(2*h); % konec