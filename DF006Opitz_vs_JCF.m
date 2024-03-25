close all
clear all
clc

format shorte

%
% Description: here I craft a Taylor interpolation at x0 of the square root
% function by sneakingly exploit JCF. Much easier, right? The only annoying part
% is I need a routine implementing th function of matrix!
%

% Define domain: I want the OPEN (a,b)
a = 0;
b = 2;
N = 100;
x = linspace( 0, 2, N + 2 ); x = x(2:N+1);

% Taylor interpolation parameters
taylor_order = 4;
x0 = 1; % center of Taylor interpolation

% operative part
sqrt_ders = @( x0,n ) ( eye( 1,n ) * sqrtm( x0 * eye( n ) + diag( ones( 1, n - 1 ),1 ) ) )'; % first row of sqrtm( J )
Tf = ( ( x(:) - x0 ).^( 0 : taylor_order-1 ) ) * sqrt_ders( x0, taylor_order );

figure,
plot( x, sqrt( x ), '-o', x, Tf(:), '-+', x0, sqrt( x0 ), '*k' )
legend('Taylor interpolant', 'Actual function', 'Approximation center','Location','SouthEast')
title('Taylor interpolation: function vs its Taylor interpolant')
axis equal


clear all % just to make sure I don't forget to instantiate anything below
%
% Description: here I craft a Newton interpolation at x0 of the square root
% function by exploiting Opitz formula.
%

% Define domain: I want the OPEN (a,b)
a = 0;
b = 2;
N = 100;
x = linspace( 0, 2, N + 2 ); x = x(2:N+1);

% Newton interpolation parameters
newton_order = 4;
z = ( cos( ( 2 * (1:newton_order) - 1 ) / ( 2 * newton_order ) * pi ) ) * ( a - b ) / 2 + ( a + b ) / 2; % interpolation points. Chebyshev in (a,b) no cap skrt skrt

% operative part
sqrt_ddfs = @( x ) ( eye( 1,numel( x ) ) * sqrtm( diag( x ) + diag( ones( 1, numel( x ) - 1 ),1 ) ) )'; % first row of sqrtm( Z )
Nf = cumprod( [ ones( size(x(:))), x(:) - z(1:end-1) ],2 ) * sqrt_ddfs( z );

figure,
plot( x, sqrt( x ), '-o', x, Nf(:), '-+', z, sqrt( z ), '*k' )
legend('Newton interpolant', 'Actual function','Interpolation points','Location','SouthEast')
title('Newton interpolation: function vs its Newton interpolant')
axis equal


disp('What happens when one or more interpolation points in Newton interpolation collapse? What are the similarities with JCF?');













%
