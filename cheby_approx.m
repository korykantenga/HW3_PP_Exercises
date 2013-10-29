function y = cheby_approx(x,theta,a,b)
	
% Chebyshev approximation of order n given a vector theta

	y = theta*chebyshev(x,0:length(theta)-1,a,b);

