close all
clear all
clc

format shorte

%
% Description: here we compare few way for computing the divided differences for
% the exponential function.
%
% You're welcome to play around with this code!
%
% Think about this: if the function to interpolate wasn't the exponential, how
% more difficult this task could have got? Would you like to play around with
% the sine function? Why the trick I used for computing the exponential's dd
% wouldn't work for the sine?
%

a = - 1;
b =   1;
n = 50;% + 100;

f = @( x ) exp( x );

interpolation_sequence_kind = 'chebyshev'; % available: random, linspaced, chebyshev, coalescing
switch interpolation_sequence_kind
  case 'random'
    z = a + rand( 1, n ) * ( b - a );
  case 'linspaced'
    z = linspace( a, b, n );
  case 'chebyshev'
    z = ( cos( ( 2 * (1:n) - 1 ) / ( 2 * n ) * pi ) ) * ( a - b ) / 2 + ( a + b ) / 2;
  case 'coalescing'
    z = [ z(3)+1e-5, z(3)+2e-5, z(3)-1e-5, z ]; % coalescing
end


'Compute reference divided differences';
d_ref = ref_dd( z );

'Standard Recurrence for divided differences';
D = diag( f( z(:) ) );
for i = 1 : n - 1
  for j = 1 : ( n - i )
    D( j, j + i ) = ( D( j, j + i - 1 ) - D( j + 1, j + i ) ) ./ ( z( j ) - z( j + i ) );
  end
end
d = D(1,:)';

'Naive Opitz application';
dd = expm( diag( z ) + diag( ones( numel( z ) - 1,1 ), 1 ) );
dd = dd(1,:)';

'Stupid Opitz application';
Z = ( diag( z ) + diag( ones( numel( z ) - 1,1 ), 1 ) );
[ V,DIAG ] = eig( Z );
Ds = V * expm( DIAG ) * inv( V );
ds = Ds(1,:)';


%
disp('Some error info');
[ d, d_ref, d - d_ref, abs( d - d_ref ) ./ abs( d ) ]
norm( d - d_ref )

max( abs( d - d_ref ) ./ abs( d ) )
max( abs( d - d_ref ) ./ abs( d_ref ) )



figure,
semilogy( 1:length(d), abs( d ), '-', 1:length(d), d_ref, '-x', 1:length(d), dd, '-o', 1:length(d), abs( ds ), '-+' )
legend('Standard Recurrence dd', 'Advanced Opitz dd', 'Naive Opitz dd', 'Stupid Opitz dd', 'Location', 'SouthWest')
%
disp('Exercise: figure out why Stupid Opitz formula is THAT stupid!')
%

%


function d = ref_dd( x )
%
% Can you work out what kind of pro-gamer move am I applying here?
%
  Z = diag( x ) + diag( ones( numel( x ) - 1,1 ), 1 );

  s = 8;
  Z = Z / 2^s;
  d = expm( Z );
  for i = 1 : s
    d = d^2;
  end
  d = d( 1,: )';

end
