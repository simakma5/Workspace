%% Úvod do prostředí MATLAB
% Revize IX 2020: Ing. Petr Krýže (kryzepet@fel.cvut.cz)
% Pro výukové potřeby předmětu B2M31AEDA a B2B31CZS, FEL, ČVUT v Praze
% -------------------------------------------------------------------------
% Tato příručka slouží pro rychlou orientaci v základech programovacího
% jazyka MATLAB. Použijte CTRL+F pro vyvolání funkce hledání, a hledejte v
% dokumentu klíčová slova, např. "vektory" nebo "tabulky".
%
% Přestože tato příručka obsahuje mnohé, používejte hlavně funkci nápovědy
% help pro jednotlivé funkce a také hledejte řešení na Googlu, kdy na
% stránách MathWorks často najdete k funkcím a úlohám dobré příklady.
% -------------------------------------------------------------------------
%% Komentáře
% Toto je komentář, protože řádek začíná znakem "%".

%{
Blokový komentář je ohraničený procentem a složenou závorkou
a může být
přes více
řádků
%}

% Dvě procenta (%%) vyznačují rozdělení na sekce. Právě zvolenou sekci
% můžete vyhodnotit (spustit) pomocí CTRL+Enter.

% Pro každou funkci je k dispozici nápověda, stačí napsat 'help nazev_prikazu'
help help

%% Práce s proměnnými
x = 1 % Přiřadí hodnotu 1 do proměnné 'x'
%  (mezery jsou volitelné, nicméně pomáhají s přehledností)

1   % Pokud k hodnotě nespecifikujete proměnnou, MATLAB uloží hodnotu do proměnné 'ans' 
    % ... "last ANSwer" - zde se ukládá i nepřiřazený výsledek výpočtů 
x   % Vypíše v konzoli hodnotu uloženou v proměnné 'x'

% Výpis hodnoty proměnné přes skript: funkce disp
disp(x)

x;  % Použití středníku zamezí výpisu příkazu do konzole - doporučeno!

% Nejdůležitější příkaz celého MATLABu! :D - clc
clc % Vyčistí vám konzoli od všeho vypsaného

save % Uloží workspace (všechny proměnné) do souboru 'matlab.mat' do aktuálního adresáře
save('nazev_souboru','x') % Uloží proměnnou 'x' do souboru 'nazev_souboru.mat'
clear x % Smaže proměnnou 'x'
clear   % Smaže všechny proměnné
load    % Načte soubor s proměnnými 'matlab.mat' z aktuálního adresáře
load('nazev_souboru') % Načte soubor s proměnnými z 'nazev_souboru.mat'

% Výpis do konzole: funkce fprintf
% Nový řádek: \n
s = 'text'; % %s je placeholder pro textové řetězce
a = 2.4869; % %.4f je placeholder pro desetinná čísla, číslo před f je 
            % počet desetinných míst k zobrazení
b = 15; % %d je placeholder pro celá čísla
fprintf("Tohle je výpis do konzole. Mohu sem vkládat %s i čísla: %.4f, %d\n", s ,a ,b)

% Pokračování v příkazu na novém řádku (pro dlouhé příkazy): ...
s = "do editoru.";
fprintf("Tohle je ukázka výpis velmi dlouhého příkazu, který se nám nevejde %s\n", ...
    s)

% Generování náhodných čísel:
M = 10; N = 10;
% 1. Náhodná čísla s normálním rozdělením (mean = 0, std = 1): 
% Funkce randn
x = randn(M,N); disp(x)
% 2. Náhodná čísla s rovnoměrným rozdělením (interval (0,1) exkluzivně):
% Funkce rand
x = rand(M,N); disp(x)
% 3. Náhodná celá čísla s rovn. rozděl. na daném intervalu (IMIN,IMAX):
% Funkce randi
IMIN = 5; IMAX = 9;
x = randi([IMIN,IMAX],M,N); disp(x)

%% Skalární aritmetika
% Operace:
% 1. Sčítání 
x = 1+2;
fprintf("Sčítání: 1+2 = %d\n", x)

% 2. Odčítání
x = 1-2;
fprintf("Odčítání: 1-2 = %d\n", x)

% 3. Násobení
x = 81*2;
fprintf("Násobení: 81*2 = %d\n", x)

% 4. Dělení
x = 13/2;
fprintf("Dělení: 13/2 = %.3f\n", x)

% 5. Mocnění
x = 2^10;
fprintf("Mocnění: 2^10 = %d\n", x)

% 6. Odmocnění
x = sqrt(1024); % Druhá odmocnina
fprintf("Odmocnění: sqrt(1024) = %d\n", x)
x = nthroot(8,3); % Libovolná - př. třetí odmocnina z 8
fprintf("Odmocnění: nthroot(8,3) = %d\n", x)

% 7. Logaritmus
x = log(1); % Přirozený logaritmus
fprintf("Logaritmus: log(1) = %d\n", x)
x = log2(64); % O základu 2
fprintf("Logaritmus: log2(64) = %d\n", x)
x = log10(1000); % O základu 10
fprintf("Logaritmus: log10(1000) = %d\n", x)

% 8. Exponenciální funkce
x = exp(1); % Eulerovo číslo
fprintf("Exponenciální funkce: exp(1) = %.3f\n", x)
x = exp(156);
fprintf("Exponenciální funkce: exp(156) = %.3e\n", x)

% 9. Zbytek po dělení (modulo)
x = mod(5,2); % Zbytek po dělení 5 číslem 2 (2*2 = 4, 5-4 = 1)
fprintf("Zbytek po dělení (modulo): mod(5,2) = %d\n", x)

%% Vektory a matice
% Definice vektorů:
% 1. Ručně vypsáním hodnot do hranatých závorek [1 2 3]
% 2. Definovat vektory můžete po sloupcích [1, 2, 3] nebo po řádcích [1; 2; 3]
% 3. Automaticky pomocí dvojtečky od:do s krokem 1 (1:10)
% 4. Automaticky definovaným krokem od:krok:do (1:0.5:10)
% 5. Funkce linspace (help linspace)

% Generování nulového vektoru: funkce zeros (help zeros)
% Generování jednotkového vektoru: funkce ones (help ones)
% Generování vektoru NaNů: funkce nan (help nan)

% Transponování vektorů (a matic)
% 1. Apostofem (') - př.: x = [1 2 3]; x = x';
% 2. Funkce transpose (help transpose)

% Adresování pozice ve vektoru/matici: 
% 1. Pro vektor x = [1 2 3 4 5];
%       Výpis třetího prvku: x(3)
% 2. Pro matici x = [1 2 3 ; 4 5 6];
%       Výpis/adresování prvku v druhém řádku a druhém sloupci: x(2,2)
% 3. Výpis/adresování posledního prvku x(end) respektive x(end,end)
% 4. Výpis/adresování všech prvků ve sloupci/řádku
%       Všechny prvky v prvním řádku: x(1,:)
%       Všechny prvky v prvním sloupci: x(:,1)

% Mazání prvku ve vektoru/matici:
% Mazání prvku na pozici 3: x(3) = []
% Mazání prvku na pozici 1,1: x(1,1) = []

% Skládání vektorů/matic:
% slozeny_vektor=[prvni_vektor druhy_vektor]
% ... pozor na rozměry vektorů (řádky vs. sloupce)!

% Operace po prvcích jsou uvozeny tečkou:
% Př.: násobení po prvcích x = [1 2 3]; y = [4 5 6]; x.*y;

% Celková suma vektoru se určí pomocí funkce sum

% Dělení matic:
% A, B jsou matice
% Zprava: x = A/B  (řešení x pro soustavu x*B = A)
% Zleva: x = A\B   (řešení x pro soustavu A*x = B)

% Determinant matice: funkce det

%% CVIČENÍ (vektory):
% Vytvořte následující vektory:
%   - Řádkový vektor 'x' obsahující hodnoty 3, 5, 1 a 7
%   - Sloupcový vektor 'x' obsahující hodnoty 3, 5, 1 a 7

% Vytvořte následující vektory bez ručního vypisování hodnot:
%   - Vektor x obsahující hodnoty od 1 do 100
%       - určete součet všech jeho prvků
%   - Řádkový vektor 'x' obsahující hodnoty od 17.3 do -3 s krokem -0.1       
%   - Řádkový vektor 'x' obsahující 100 nul
%   - Sloupcový vektor 'x' obsahující 100 jedniček
%   - Vektor 'x' obsahující na pozici 2 až 99 samé jedničky, na pozici 1 a 100 jsou nuly

% Vytvořte vektor 'v1' obsahující hodnoty 1, 2 a 3 a vektor 'v2' obsahující
% hodnoty 2, 3 a 4. Proveďte následující operace:
%   - K vektoru 'v1' přičtěte číslo 1 a poté ho vynásobte dvěmi
%   - Vynásobte mezi sebou po prvcích 'v1' a 'v2'
%   - Umocněte po prvcích 'v2' na třetí

%% CVIČENÍ (matice):
% Vytvořte následující matici:
% A = [ 1 2 3 
%       4 5 6 
%       7 8 9 ];

% - Vypište hodnotu prvku nacházejícího se uprostřed matice 'A' (pozice 2,2) 
% - Změnte hodnotu tohoto prvku na 0
% - Vypište všechny hodnoty ve druhém řádku matice 'A'
% - Vypište všechny hodnoty ve třetím sloupci matice 'A'

% Vytvorte řádkový vektor 'x' obsahující hodnoty 3, 5 a 1
% a řádkový vektor 'y' obsahující hodnoty 2, 8 a 4.
%   - Proveďte maticové násobení mezi 'x' a transponovanou-'y'
%   - Vytvořte matici A spojením 'x' a 'y' po řádcích (dejte je na sebe)
%   - Vytvořte matici A spojením 'x' a 'y' po sloupcích (dejte je vedle
%   sebe)

% Vytvořte následující matici:
% E = [ 2 3 4 
%       4 3 2 
%       4 4 1 
%       1 2 1 ];
%   - Vytvořte matici F obsahující všechny hodnoty z E vyjma třetího řádku
%   - Vytvořte matici G obsahující všechny hodnoty z E vyjma druhého řádku
%   - Přičtěte matici F k matici G
%   - Pronásobte matici F a matici G

% Mějme následující soustavu rovnic:
%     |3  -3   4|   |U1|   |30| 
%     |1   6   5| . |U2| = | 7| 
%     |1  -2   3|   |U3|   |17|
% Určete U1, U2, U3.
% (Help: https://cs.wikipedia.org/wiki/Cramerovo_pravidlo)

%% Tabulky
% Tabulky umožňují uchovávat proměnné různých formátů (v sloupcích) pro
% příslušná pozorování (řádky). Tabulku lze vytvorit pomocí funkce table:
A = {'Petr' ; 'Pavel' ; 'Marie'};
B = [1 ; 2 ; 3];
C = [2 ; 3 ; 1];
D = [3 ; 1 ; 2];

T = table(A,B,C,D, ...
    'VariableNames',{'Jmeno','Matematika','Dejepis','Krasopis'});
disp(T)

% Někdy je praktické jednotlivé řádky "pojmenovat":
T = table(B,C,D,'RowNames',A,...
    'VariableNames',{'Matematika','Dejepis','Krasopis'});
disp(T)

% Převedení tabulky na matici:
% Je nutné volit sloupce se stejným datovým typem
x = table2array(T(:,2:end));
disp(x)
x = T{:,2:end}; % alternativně
disp(x)

% Vytvoření tabulky z matice:
T1 = table([B C D]); % Pozor! Celá matice je v rámci proměnné
disp(T1)
T2 = array2table([B C D]); % array2table vytvoří tabulku po sloupcích matice
disp(T2)

% Skládání tabulek:
T = [table(A) T2];
disp(T)

% Přístup a změna sloupců (variables):
T.Properties.VariableNames={'Jmeno','Matematika','Dejepis','Krasopis'};
disp(T)

% Přístup a změna dat v tabulce:
disp(T.Dejepis); % Sloupec dějepis
disp(T([1 2],3)); % Řádek 1 a 2 ze 3. sloupce
T.Dejepis = [1 ; 5 ; 4]; % Přiřadí hodnoty 1, 5 a 4 jako známky z dějepisu

%% CVIČENÍ (tabulky):
% Z tabulky T vymažte Dějepis a za Krasopis přidejte AED. Známky AED
% vypočtěte automaticky, aby průměr všech známek každého žáka byl 2.

%% Logické funkce
x = false; y=true; % Hodnoty log. 0 a 1

% 1. Negace (NOT)
z = ~x; disp(z);

% 2. Logický součin (AND)
z = x & y; disp(z);
z = and(x,y); disp(z); % alternativně

% 3. Logický součet (OR)
z = x | y; disp(z);
z = or(x,y); disp(z); % alternativně

% 4. Exkluzivní logický součet (XOR)
z = xor(x,y); disp(z);

% 5. Any - vrací TRUE je-li alespoň jeden prvek logická 1 (TRUE)
x = [false true false false];
z = any(x);  disp(z);

% 6. Porovnávání
a = [1 2 3 4 5 6 7 8 9 10];
% Vrací logické hodnoty true (1) a false (0) podle výsledku srovnání
z = a < 3; disp(z); % Menší
z = a >= 3; disp(z); % Větší rovno
z = a <= 3; disp(z); % Menší rovno

z = a == 5; disp(z); % Vrátí log. 1 pro prvky rovné 5
z = a ~= 5; disp(z); % Vrátí log. 1 pro prvky nerovné 5

% 7. Hledání
a = [1 5 3 4 5 6 7 8 9 5];
z = find(a == 5); % Vypíše polohy (indexy) prvků rovných 5

a = [1 5 3 ; 4 5 6 ; 7 8 5];
[row, col] = find(a > 5); % Vypíše souřadnice prvků s hodnotou větší než 5

%% Podmínky IF, ELSE, SWITCH
% 1. Napíše 'AED cant melt steel beams' vždy (triviální):
if 1 % Alternativně lze napsat 'if true'
    disp('AED cant melt steel beams')
end 

% 2. Porovná hodnoty a podle výsledku vypíše:
podminka = 1;
if (podminka ~= 0)
    disp('AED cant melt steel beams')
end 

% 3. ELSE
podminka = false;
if podminka
    disp('This statement is false.'); 
else
    disp('This statement is true.');
end 

% 4. SWITCH
jmeno = 'petr';
switch jmeno
    case 'martin'
        fprintf("Jméno je Martin.\n")
    case 'petr'
        fprintf("Jméno je Petr.\n")
    case 'jan'
        fprintf("Jméno je Jan.\n")
    otherwise
        fprintf("Jméno je neznámé...\n")
end

%% CVIČENÍ (logické funkce a podmínky):
% Vygenerujte šachovnici 'R' o rozměrech 1000x1000, kde logická 0 bude
% reprezentovat tmavá pole a logická 1 bude reprezentovat světlá pole.
% Pokuste se nepoužívat žádné cykly.

% Napište následující podmínku:
% Jestli je v předchozí matici 'R' na první pozici 1, potom všechny prvky 
% matice logicky invertujte. V opačném případě vypište "Trefil jsem to dobře!"

%% Cykly FOR, WHILE
% 1. FOR cyklus
for i = 1:10 % Cykluje od jedničky do 10 (včetně), hodnotu přiřadí do 'i'
    disp(i)
end

k = 1:10; % Procházení předem existujícího vektoru 
for i = 1:length(k) % Length - počet prvků v 1-D vektoru/matici
    disp(k(i))
end

% Vnořené for cykly (procházení 2D prostoru)
m = 5; % Počet řádků
n = 5; % Počet sloupců
X = randn(m,n); % Náhodná matice
for i = 1:m % Prochází řádky
    for j = 1:n % Prochází sloupce
        fprintf("Řádek: %d, Sloupec: %d, Hodnota: %.3f\n", i, j, X(i,j))
    end
end

% 2. WHILE cyklus
i = 0; 
while (i < 11) % Dokud bude platit podmínka v argumentu, cyklus bude probíhat
    i = i + 1;
    disp(i);
end

%% CVIČENÍ (for a while cykly):
% Zkuste si vygenerovat šachovnici (jako v předchozí části) s pomocí for cyklů.

% Zkuste si vygenerovat šachovnici znovu, tentokrát pomocí while cyklu.

%% Funkce
% Lze je napsat do úplně samostatného souboru (který pak musí být součástí
% MATLABovské cesty (PATH), aby byl spustitelný. Soubor se pak musí
% jmenovat stejně jako název funkce. (funkce.m)

% Také lze napsat si dílčí funkce přímo do skriptu, na jeho konec (nutné).
% Funkce lze pak libovolně ve skriptu používat (což je fajn pro mnohokrát
% opakované části kódu)

% 1. Definice:
% Př.: Funkce na manuální násobení
% x je výstupní hodnota (té přidělíme v těle funkce hodnotu - ekvivalent return
% 'a' a 'b' jsou vstupy, může jich být libovolný počet
% PODÍVEJTE SE NA KONEC TOHOTO SKRIPTU

% 2. Použití funkce - triviální:
x = vynasob(5,6); disp(x);

% 3. Kontroly počtu argumentů:
% Funkce nargin - vrátí počet argumentů které byly předány na vstupu funkce
% Funkce narginchk - zkontroluje správný počet argumentů (podrobněji viz
% help)

%% CVIČENÍ (funkce):
% Zkuste přepsat kód generování šachovnice do jedné univerzální funkce...
% chessboard = generate_chessboard(N)
% ...kde chessboard je výsledná šachovnice velikosti NxN a N je délka
% strany čtvercové šachovnice

%% Příklad funkce
function x = vynasob(a,b)
    x = a*b;
end



