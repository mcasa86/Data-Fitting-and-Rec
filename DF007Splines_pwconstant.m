close all
clear all
clc
format shorte

%
% Description: here we play around with piecewise constant interpolation of some
% interpolation data ( x_i, y_i ).
%
% You're welcome to play around with this code!
%


N = 1e4; % we want to avoid stupid artifacts in PLOTTING
n = 50;  % number of knots, i.e., interval breaks

a = - 1;
b =   1;
x     = linspace( a,b, N );
knots = linspace( a,b,n ); % sort( rand( 1,n ) * ( b - a ) + a );
knots(   1 ) = knots(   1 ) - eps; % could you say why did I do this?
knots( end ) = knots( end ) + eps; % could you say why did I do this?

p = 1 / 5; % any convex combination is fine
% btw one usually chooses knots AFTER the interpolation points
interp_points = p * knots( 1 : end-1 ) + ( 1 - p ) * knots( 2 : end );
y = cos( interp_points * 2 * pi ); %randn( 1,length( knots ) );

figure,
for i = 2 : n
  ids = find( ( x >= knots( i - 1 ) ) .* ( x <= knots( i     ) ) );
  plot( x( ids ), y( i - 1 ) * ones( size( ids ) ), '--k' ), hold on
end
xlim([ a,b ])
ylim([ min( y ),max( y ) ])
plot( x, 0 * x, ':k'), hold on
plot( knots, 0 * knots, 'xk', 'MarkerSize', 5, 'Linewidth', 2 ), hold on
plot( interp_points, 1 *     y, 'ok', 'MarkerSize', 5, 'Linewidth', 2 ), hold on


%
con_spline = zeros( size( x ) );
for i = 2 : n
  ids = find( ( x >= knots( i - 1 ) ) .* ( x <= knots( i     ) ) );
  con_spline( ids ) = y( i - 1 ) * ones( size( ids ) );
end
% figure,
ids = find( ( x >= knots( 1 ) ) .* ( x <= knots( end ) ) );
plot( x( ids ), con_spline( ids ), '-', 'Linewidth', 2 ), hold on
xlim([ a,b ])
ylim([ min( y ),max( y ) ])



%
% Splines: motivation in the context of data fitting and reconstruction. Least squares refresher, average as constant least square best approximant of a set of values. Piecewise constant splines; the linear space $1. Uniqueness of the piecewise constant interpolant spline. Upper bound for uniform error committed using piecewise constant splines (with proof). Mean value theorem for divided differences. Fourier decomposition refresher and light introduction to the topic of the first project: the Haar Wavelets.
