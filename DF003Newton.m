close all
clear all
clc

format shorte

%
% Description: here I craft a Newton interpolation at a interpolation sequence
% of the sine function. How different is it from Taylor interpolation?
%
% You're welcome to play around with this code!
%
% Think about this: if the function to interpolate wasn't the exponential, how
% more difficult this task could have got?
% If the interpolation point coalesced, how quick accuracy would go down the
% toilet?
%

a = -2 * pi;
b =  2 * pi;
n = 5;
N = 100;

f = @( x ) exp( x );

interpolation_sequence_kind = 'chebyshev'; % available: random, linspaced, chebyshev
switch interpolation_sequence_kind
  case 'random'
    z = a + rand( 1, n ) * ( b - a );
  case 'linspaced'
    z = linspace( a, b, n );
  case 'chebyshev'
    z = ( cos( ( 2 * (1:n) - 1 ) / ( 2 * n ) * pi ) ) * ( a - b ) / 2 + ( a + b ) / 2;
end

% evaluation points
x = linspace( a, b, N );

% compute divided differences
D = diag( f( z(:) ) );
for i = 1 : n - 1
  for j = 1 : ( n - i )
    D( j, j + i ) = ( D( j, j + i - 1 ) - D( j + 1, j + i ) ) ./ ( z( j ) - z( j + i ) );
  end
end

Nf = zeros( size( x(:) ) ); % this is going to be the evaluation of T[f](x) at the evaluation points
acc = ones( size( x(:) ) ); % we need this for performance
for i = 1 : n
  Nf = Nf + D( 1,i ) * acc;
  acc = acc .* ( x(:) - z( i ) ); % here I update it once too much, not a big problem come on
end


figure,
plot( x, Nf, 'o', x, f( x ), '-' )
legend('Newton interpolant', 'Actual function')
title('Newton interpolation: function vs its Newton interpolant')



%
