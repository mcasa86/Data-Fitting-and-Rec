close all
clear all
clc

format shorte

% choose your favourite eigenvalues
lambda = 1:4;

% generate a matrix with controlled eigenvalues
n = length( lambda );
A = randn( n );
A = A * A'; % symmetric matrix = real eigenvalues and eigenvectors
[ V,D ] = eig( A ); % generated n distrinct eigenvectors
D = diag( lambda ); % replace eigenvalues with eigenvalues of our choice :3

A = V * D * inv( V ); % A = V * D * V';

if any( not( diff( sort( diag( D ) ) ) ) )
  warning('Make sure you chose simple eigenvalues or everything works fine apart for the sine of matrix down low.');
end

'the funniest thing';
%
% I want to produce two DIFFERENT polynomial p & q of degree N > n that take the
% same values at lambda!
%
% What values? Say sin( lambda ), why not? But, to be different, p & q have to
% assume different values at the subset of N - n values z.
% Say one takes values sin( z ) and the other exp( z ) / 20 + 1 so they're
% different for sure, okay?
%
N = n + 3; % don't exaggerate
z = randn( 1, N - n );

p_at_interpolation_sequence = [ sin( lambda ), sin( z )          ];
q_at_interpolation_sequence = [ sin( lambda ), exp( z ) / 20 + 1 ];

% let's be brave and use the infamous Vandermonde interpolation technique
interpolation_sequence = [ lambda, z ];

Vandermonde = fliplr( vander( interpolation_sequence ) ); % could you explain why the vandermonde matrix is often ill conditioned?

coefficients_of_p = Vandermonde \ p_at_interpolation_sequence(:);
coefficients_of_q = Vandermonde \ q_at_interpolation_sequence(:);

figure,
x = linspace( min( interpolation_sequence ), max( interpolation_sequence ) );
plot( x, polyval( flip(coefficients_of_p),x ), '-', x, polyval( flip(coefficients_of_q),x ), '--' )
hold on
plot( lambda, sin( lambda ), '*k' )
axis equal

p_of_A   = polynomial_evaluation( A, coefficients_of_p );
q_of_A   = polynomial_evaluation( A, coefficients_of_q );
sin_of_A = my_sinm( A );

disp('p_of_A');
p_of_A
disp('q_of_A');
q_of_A
disp('sin_of_A');
sin_of_A

function p_of_X = polynomial_evaluation( X, coefficients )
  p_of_X = coefficients( 1 ) * eye( size( X ) );
  for i = 1 : length( coefficients ) - 1
    p_of_X = p_of_X + coefficients( i + 1 ) * X^i; % not efficient but we don't care, do we?
  end
end

function sA = my_sinm( A )
  [ V,D ] = eig( A );
  sA = V * diag( sin( diag( D ) ) ) * inv( V ); % A is meant to have simple eigenvalues, why didnt i simply use jordan canonical form?
end



%
% Some references I find interesting:
%
% Yousef Saad's "Iterative Methods for Sparse Linear Systems", 2003.
% Nicholas J. Higham's "Functions of Matrices - Theory and Computation", 2008.
% Carl de Boor, "Divided differences". Surveys in Approximation Theory, 1:46–69, 2005.
