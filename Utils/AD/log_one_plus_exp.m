function y = log_one_plus_exp(x,delta)
%An AD friendly version of log(1+delta*exp(x))
%Details from Mark Greenstreet April 12th 2018
%
%log(x) = u + log(x*exp(-u) for any u
%therefore:
%log(1+x) = x + log((1+x)*exp(-x)), where u = x
%
%log(1+delta*exp(y)) = exp(y) + log((1+exp(y))*exp(-exp(y))) 
%with y = exp(x)
%
%log(1+delta*exp(x)) =
%delta*exp(x)+log((1+delta*exp(x))*exp(-exp(delta*x)))
%
%log(1+exp(x)) = x + log((1+exp(x))*exp(-x) = x+log(1+exp(-x))

y = 0*x; % allocate an AD value
i_neg = find(x<=0);
i_pos = find(x>0);
y(i_neg) = exp(x(i_neg)) + log(exp(-exp(x(i_neg)))+exp(x(i_neg)).*exp(-exp(x(i_neg))));
y(i_pos) = x(i_pos) + log(1.0 + (delta + 1.0) .* exp(-x(i_pos)));
end

