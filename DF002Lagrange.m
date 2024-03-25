close all
clear all
clc

format shorte

%
% Description: here I astutely craft a Lagrange interpolation at a interpolation
% sequence of the sine function. How different is it from Taylor interpolation?
%
% You're welcome to play around with this code!
%
% Think about this: if the function to interpolate wasn't the exponential, how
% more difficult this task could have got?
%

a = -2 * pi; % left extreme
b =  2 * pi; % right extreme
n = 5;   % number of interpolation points
N = 100; % number of evaluation points

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

% can you spot the funny trick I've used?
Z = z(:) - z(:)';
X = x(:) - z(:)';
Lf = zeros( size( x(:) ) ); % this is going to be the evaluation of T[f](x) at the evaluation points
for i = 1 : n
  ell{ i } = prod( X( :,[ 1:i-1, i+1:end ] ) ./ Z( i,[ 1:i-1, i+1:end ] ), 2 );
  Lf = Lf + ell{ i } .* f( z( i ) );
end

figure,
plot( x, Lf, 'o', x, f( x ), '-' )
legend('Lagrange interpolant', 'Actual function')
title('Lagrange interpolation: function vs its Lagrange interpolant')

disp('Pointwise error:');
[ Lf - f( x(:) ) ]

figure,
plot( x, ell{ 3 }, '-')
title('Plot of \ell_3(x)')
axis equal

figure,
l = 0;
for i = 1 : length( ell )
  l = l + ell{ i };
end
plot( x, l, '-')
title('Plot of the sum of all \ell_i(x), i = 0,1...')
axis equal




%
