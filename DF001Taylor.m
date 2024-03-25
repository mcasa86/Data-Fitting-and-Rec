close all
clear all
clc

format shorte

%
% Description: here I astutely craft a Taylor interpolation at x0 of the sine
% function.
%
% You're welcome to play around with this code!
%
% Think about this: if the function to interpolate wasn't the sine, how more
% difficult this task could have got?
%

% Matlab has this funny data structure called "cell". Run "help cell" in the
% Command Window to have more information. This is very useful to store data
% that is heterogeneous or, roughly speaking, not numbers.
% In f{ 1 } we store the function we want to approximate, i.e., sin( x ).
% In f{ i } we store its (i-1)th derivative. We need this to call the right
% (cyclic) derivative right when we need it!
f{ 0 + 1 } = @( x )   sin( x );
f{ 1 + 1 } = @( x )   cos( x );
f{ 2 + 1 } = @( x ) - sin( x );
f{ 3 + 1 } = @( x ) - cos( x );


taylor_order = 30;
x0 = 5; % center of Taylor interpolation
a = - 2 * pi; % left extreme
b =   2 * pi; % right extreme


N = 100+1; % arbitrarily chosen number of evaluation points
x = linspace( a,b, N )'; % generate N points from our interval

Tf = zeros( size( x ) ); % this is going to be the evaluation of T[f](x) at the evaluation points
acc = ones( size( x ) ); % we need this for performance
for i = 0 : taylor_order
  cyclic_id = mod( i, 4 ) + 1; % try and make sense of this!!
  Tf = Tf + f{ cyclic_id }( x0 ) / factorial( i ) .* acc;
  acc = acc .* ( x - x0 );
end

figure,
plot( x, f{ 0 + 1 }( x ), '-', x, Tf, '-+', x0, f{ 0 + 1 }( x0 ), '*k' )
legend('Taylor interpolant', 'Actual function')
title('Taylor interpolation: function vs its Taylor interpolant')
axis equal


disp('Pointwise error:');
[ Tf - f{ 0+1 }( x ) ]

figure,
semilogy( x, abs( [ Tf - f{ 0+1 }( x ) ] ./ ( f{ 0+1 }( x ) + realmin ) ) )
hold on
plot( x0, 1e-15, 'k*')
title('Logarithmic plot of the error with Taylor center highlighted, what do you notice?')
