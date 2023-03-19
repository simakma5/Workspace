%Function CDE: Euclidan cepstral distance between two cepstra.
%
%Usage: D=cde(c1,c2,cp);
%
%  Mandatory input parameters:
%    c1  - vector/matrix of cepstral coefficients of the first spectrum
%    c2  - vector/matrix of cepstral coefficients of the second spectrum
%          (particular cepstra must be in rows of matrices c1 and c2)
%
%  Optional input parameters:
%    cp  - number of used cepstral coefficients     <minimum of size(c1,2),size(c2,2)>
%
%  Output: D  - vector (value) of Euclidan cepstral distance 
%
%    1) when c1 and c2 are row vectors of same dimensions, D is Euclidan distance 
%       between corresponding cepstra in rows of matrices c1 and c2
%
%    2) when c1 and c2 are matrices of same dimensions, D is Euclidan distance 
%       between corresponding cepstra saved in rows of matrices c1 and c2
%
%    3) when c1 or c2 is a vector, D is Euclidan distance this vector and
%       all cepstra in rows of the second matrix 
%
%    4) when the number of columns (i.e. a number of particular cepstral coefficients)
%       of c1 and c2 is not same, the distance is computed on the basis of minimum
%       available cepstral coefficients

%                           Made by Petr Pollak
%                              26 Oct 2022
%                         CVUT FEL K13131, Prague

function D=cde(c1,c2,cp);

if nargin<2,
  error('ERROR in function CDE: Not enough input arguments!');
end

size_c1=size(c1);
size_c2=size(c2);

if nargin<3,
  cp=min([size_c1(2) size_c2(2)]);
  if ( size_c1(1) ~= size_c2(1) ),
    disp('Warning in function CDE: Different numbers of cepstral coefficients in c1 and c2.');
    disp('  Number of cepstral coefficients was set according to smaller value!');
end
end;

if cp>min([size_c1(2) size_c2(2)]);,
  disp('ERROR in function CDE: Inconsistent number of cepstral coefficients.');
  error('  Check dimensions of c1 or c2!');
end

c1=c1(:,1:cp);
c2=c2(:,1:cp);


if ( size_c1(1) == 1 | size_c2(1) == 1 ),

    if ( size_c1(1) == 1 ),
        c0=c1;
        cc=c2;
    else
        c0=c2;
        cc=c1;
    end
    D=sqrt(sum((cc-c0).^2,2));;
    
    return

elseif ( size_c1(1) == size_c2(1) ),

    D=sqrt(sum((c1-c2).^2,2));

    return  

else

    error('Inconsistent number of cepstral coefficients.\n  Check dimensions of c1 or %s!','c2');

end


