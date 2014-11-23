function mm = quadConvert(m,src,dst)
% mm = quadConvert(m,src,dst)
%  The function converts between quadratic matrix
%  'b': the coefficient of quadratic polynomial as a vector
%  'B': the matrix form of 'b', which is a lower triangular matrix
%  'Bfull': copy lower triangular to upper triangular except the diagnoal elements
%  'A": the symmetric matrix form of quadratic polynomial
% 
% m: each column is one matrix to convert 
% src,dst: matrix form
% mm: the result, column for a cube.

[k,nc] = size(m);

% conver to A matrix first
switch(lower(src))
case lower('b')
  b = m;
  n = round((sqrt(8*k+1)-3)/2);
  ind = tril(true(n+1,n+1));
  map = zeros(n+1,n+1);
  map(ind) = 1:k;
  map = map+map'-diag(diag(map));
  Bfull = b(map,:);  
  A = repmat(reshape((1+eye(n+1))/2,[],1),1,nc).*Bfull;
case lower('B')
  B = m;
  n = round(sqrt(k)-1);
  mB = reshape(B,n+1,[]);
  tB = reshap(mB',n+1,[]);
  Bfull = (mB+tB)/2;
  Bfull = reshape(Bfull,[],nc);
  A = repmat(reshape((1+eye(n+1))/2,[],1),1,nc).*Bfull;
case lower('Bfull')
  Bfull = m;
  n = round(sqrt(k)-1);
  A = repmat(reshape((1+eye(n+1))/2,[],1),1,nc).*Bfull;
case lower('A')
  A = m;
  n = round(sqrt(k)-1);
otherwise
  error('do not support now');
end

% convert A to dst
switch(lower(dst))
case lower('b')
  Bfull = reshape(reshape((2-eye(n+1)),[],1),1,nc).*A;
  ind = tril(true(n+1,n+1));
  b = Bfull(ind,:);
  mm = b;
case lower('B')
  Bfull = reshape(reshape((2-eye(n+1)),[],1),1,nc).*A;
  B = zeros(k,nc);
  ind = tril(true(n+1,n+1));
  B(ind,:) = Bfull(ind,:);
  mm = B;
case lower('Bfull')
  Bfull = reshape(reshape((2-eye(n+1)),[],1),1,nc).*A;
  mm = Bfull;
case lower('A')
  % done
  mm  = A;
otherwise 
  error('do not support now');
end;
