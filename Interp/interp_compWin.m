function c = interp_compWin(n)
% Compute cosine windows for mspice that makes the interplation c^n continuous.
%   n should be a odd number. (same result for n=2*m and n=2*m+1)
% A cosine window is of the form 
%   w(x) = sum_{i=0}^{n} c_i * cos(i*x*pi)
%
% The first constraint is that the sum of weight is 1:
%   I. w(x)+w(1-x) = 1; for all x \in [-1,1]
% Thus
%   w(x)+w(1-x) = 2( c_0 + sum_{i=2:2:n} c_i * cos(i*x*pi) 
%   => c_0 = 1/2, c_{2k} = 0, for all k \in [1:2:n/2]
% Therefore, the consine window is 
%   w(x) = 1/2 + sum_{i=1:2:n} c_i * cos(i*x*pi)
% ( n must be an odd number as c_{2*i} = 0. Let n=2*m+1)
% 
% Other constraints include 
%   II. w(0) = 1  (weight is 1 at the center)
%   III. (d^k)/(dx^k) w(1) = 0, for all k \in [0,n]. (c^d continuous)
% II => sum_{i=1:2:n} c_i = 1/2 (eq1) 
% III.  When k is odd => (d^k(/(dx^k) w(1) = ... * sin(i*x*pi) = 0
%   Therefore, only care about k = 2:2:n  ( w(1)=0 when k==0, same with eq1) 
%     (d^k)/(dx^k) w(x) =  
%       sum_{i=1:2:n} c_i * i^k * (-1)^(k/2) * cos(i*x*pi)
%   =>  sum_{i=1:2:n} c_i * i^k * (-1)^(k/2) * -1 = 0, for all k = 2(1:m)  (eq2)
% Therefore, we have m+1 variables and m+1 conditions. 
% The value of a can be solved by a linear system (eq1,eq2)

m = floor(n/2);
A = zeros(m+1,m+1);
for k=2*(0:m)
  row = k/2+1;
  for i=1:2:2*m+1
    col = (i+1)/2;
    A(row,col) = i^k * (-1)^(k/2);
  end
end
b = [1/2;zeros(m,1)]; % first element for w(1)=0 (eq1)
c = A\b;
s = sprintf('The cosine window is \n 0.5');
for i=1:m+1
  if(c(i)>0), op = '+'; else op = ''; end
  s = sprintf('%s%s%f*cos(%d*x*pi)',s,op,c(i),2*i-1);
end
disp(s);
c = [0.5;c];
