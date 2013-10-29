function y = col_points(a,b,n)

% Collocation points for chebyshev polynomial of order n between a and b

	y = cos((2*(n:-1:1)-1)*pi/(2*n));
	y = (y+1)*(b-a)/2+a;

