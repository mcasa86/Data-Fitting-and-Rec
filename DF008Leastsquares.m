close all
clear all
clc

format shorte


x = linspace( 0,1 );
noise = randn( size( x ) ) / 10;

'EXAMPLE 01';
y = @( x ) 3 * x + 2;

noisy_data = y(x(:)) + noise(:);

A = [ ones( size( x(:) ) ), x(:) ]; % this is a vandermonde at the end of the day

coefficients_babies = A \ noisy_data; % coefficients_dumbies  = inv( A ) * noisy_data;
coefficients = ( A' * A ) \ ( A' * noisy_data );

figure,
plot( x, y(x) + noise, '*', x, y(x), '.', x, coefficients( 1 ) + coefficients( 2 ) * x,'-' )
legend('noisy','actual','fitted','Location','SouthEast')

% return

'EXAMPLE 02';
y = @( x ) - 3 * x.^2 + x - 1;

A = [ ones( size( x(:) ) ), x(:), x(:).^2 ]; % this is a vandermonde at the end of the day

noisy_data = y(x(:)) + noise(:);
coefficients = ( A' * A ) \ ( A' * noisy_data );


figure,
plot( x, y(x) + noise, '*', x, y(x), '.', x, coefficients( 1 ) + coefficients( 2 ) * x + coefficients( 3 ) * x.^2,'-' )









%
