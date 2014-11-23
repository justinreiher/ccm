function [b, ibnds, A, u,vbnds] = lft_simu_fit(M, v)
% [b, ibnds] = lft_simu_fit(M, v)
% This function works for any n-d matrix M and n-d cube v.
%
%  v specifies the domain cube for which we are producing a fit:
%    v(i,1) is the lower bound for domain variable i.
%    v(i,2) is the upper bound for domain variable i.
%  Let m = size(v,1)+1.  The extra one is for a "variable" that we
%    add which is 1 everywhere -- this provides the constant term for
%    the regression.
%  M.YC is an array with
%    M.YC(i1, i1, ... im) =
%        SUM_{j1=1}^i1 SUM_{j2=1}^i2 ... SUM_{jm=1}^im y(j1, j2, ... jm)
%  M.GRID.v0 and M.GRID.dv are vectors such that the value of the i^th element for
%    domain variable j is M.GRID.v0(j) + j*M.GRID.dv(j).
%
% Return values:
%  b: the best least-squares fit in the cube specified by v.
%       yy = b'*[vv; 1]
%     where yy is the estimate for the current, and vv is the
%     voltage vector where the estimate is desired.  length(vv) = size(v,1).
%     The last element of b is the constant term from the regression,
%     hence the appending of a 1 to vv.
%  ibnds: indices describing the cube in M.YC used for the regression.
%     The j^th dimension of this cube ranges from ibnds(j,1) to ibnds(j,2).

% Derivation of the least-squares optimization:
% min_b SUM_{k=1}^n (y(k) - yy(k))^2
%   Where yy(k) = b^T x(k), and b and x are vectors of length m,
%   and n is the number of points that we are fitting.
% Let z = SUM_{k=1}^n ||y(k) - yy(k)||^2
% Taking derivatives, we get:
%   dz/db_i = d/db_i (SUM_{k=1}^n (y(k) - SUM_{j=1}^m b_j * x_j(k))^2)
%      = SUM_{k=1}^n ( 2 * ((SUM_{j=1}^m b_j * x_k(k)) - y(k))
%            * d/db_i (Sum_{j=1}^m (b_j * x_j(k)))
%                       )
%      = SUM_{k=1}^n 2 * x_i(k) * ((SUM_{j=1}^m b_j * x_j(k)) - y(k))
%      =   (2 * SUM_{j=1}^m b_j * SUM_{k=1}^n (x_i(k) * x_j(k)))
%        - (2 * SUM_{k=1}^n x_i(k) * y(k))
% Let A be the matrix with A(i,j) = (SUM_{k=1}^n (x_i(k) * x_j(k))))/n.
% Let u be the vector with u(i) = (SUM_{k=1}^n x_i(k) * y(k))/n
% z is minimized when dz/db = 0.  This means:
%   b = A \ u

% First, we find the indices for the cube.
% We take floor or ceil to make sure that our subcube of the sample
% contains the cube described by v (the components of v are continuous).
% ibnds(i,1) is the lowest index for the i^th dimension of M.YC for our cube.
% ibnds(i,2) is the highest index for the i^th dimension of M.YC.
% vbnds is the approximation of v corresponding to ibnds.
if(nargin<2 || isempty(v)) % use all, maximum bbox.
  nv = size(M.GRID.v0,1);
  v0 = [M.GRID.v0, M.GRID.v0]; dv = [M.GRID.dv, M.GRID.dv];
  ibnds = ones(nv,2);
  ibnds(:,2) = size(M.YC)';
  vbnds = [v0+dv.*ibnds;[1,1]];
else
  nv = size(v, 1); %assert(nv==size(M.GRID.v0,1));
  v0 = [M.GRID.v0, M.GRID.v0];  dv = [M.GRID.dv, M.GRID.dv];
  ibnds = (v - v0) ./ dv;
  ibnds = [floor(ibnds(:,1)), ceil(ibnds(:,2))];
  % fix bug, exceed boundary of model because of round off
  ibnds = [max(1,ibnds(:,1)), min(M.GRID.nv,ibnds(:,2))];
  vbnds = [ v0 + dv .* ibnds; [1,1]];
end;

% Calculate the u vector.
% In the code below, lft_simu_fit_ycv_help(M.YC, ibnds, i, j, nv, 1) returns a vector
% that goes along the % j^th dimension of M.YC.
% Let
%    yj(i) =
%      SUM_{i1=lo(1)}^hi(1) SUM_{i2=lo(2)}^hi(2) ... SUM_{i_nv=lo(nv)}^hi(nv)
%   y(i1, i2, ... i, ... im)
% where the sums skip dimension j, and the i in y(i1, i2, ... i, ... im) is
% for dimension j.  lo(i) and hi(i) are the lower and upper bounds for
% indices for dimension i; in other words, lo(i) = ibnds(i,1), and
% hi(i) = ibnds(i,2).  Likewise, v_i(lo(i)) = vbnds(i,1), and
% v_i(hi(i)) = vbnds(i,2).  Now,
%   ycv(i) = sum_{k=1}^{i-1} yj(lo(i) + k - 1)
% with ycv(i) defined for i=1..(ibnds(j,2) - ibnds(j,1) + 2).
% We want to calculate sum_{k=lo(j)}^hi(j) vj(k) * yj(k).
% Note that yj(lo(j)+k) = ycv(k+1) - ycv(k), and vj(k) = M.GRID.v0(j) + k*M.GRID.dv(j).
% Let n(j) = hi(j) - lo(j) + 1.
% We get:
%     sum_{k=lo(j)}^hi(j) vj(k) * yj(k)
%   = sum_{k=1}^n(j) vj(lo(j) + k - 1) * (ycv(k+1) - ycv(k))
%   =   vj(hi(j))*ycv(n(j)+1)
%     + (sum_{k=2}^n(j) ycv(k)(vj(lo(j)+k-2) - vj(lo(j)+k-1)))
%     - vj(lo(j))*ycv(lo(j))
%   =   vbnds(j,2)*ycv(hi(j)+1)
%     - M.GRID.dv(j) * (sum_{k=2}^n(j) ycv(k))
%     - vbnds(j,1)*ycv(1)
% We divide this by prod(n) to get u(1:nv).
% Finally, we need u(nv+1).  This is just the sum of all of the y's in the
% cube.  We can calculate this from any ncv vector:
%    u(nv+1) = ycv(end) - ycv(1).
%  There's one detail left.  The last row of vbnds is for the
% constant term.  It's not represented in ibnds because M.YC doesn't
% have a corresponding dimension.  We treat it as a singleton dimension
% by appending a 1 to the n vector.
u = zeros(nv+1, 1);
i = cell(nv, 1);
for j = 1:nv
  % replace 
  %i{j} = (ibnds(j, 1)-(ibnds(j,1) > 1)):(ibnds(j,2));
  %ycv = lft_simu_fit_ycv_help(M.YC, ibnds, i, j, nv, 1);
  % with
  if(ibnds(j,1)>1)
    i{j} = (ibnds(j,1)-1) : ibnds(j,2);
    ycv = lft_simu_fit_ycv_help(M.YC,ibnds,i,j,nv,1);
  else
    i{j} = ibnds(j,1):ibnds(j,2);
    ycv = lft_simu_fit_ycv_help(M.YC,ibnds,i,j,nv,1);
    ycv = [0;ycv]; % insert 0 for ycv(ibnds(j,1)-1);
  end;
  assert(length(ycv)==ibnds(j,2)-ibnds(j,1)+1+1);% nj+1
  u(j) = vbnds(j,2)*ycv(end) - M.GRID.dv(j)*sum(ycv(2:end-1)) - vbnds(j,1)*ycv(1);
  if(j==nv)
    u(nv+1) = ycv(end) - ycv(1);
  end;
end;
n = [(diff(ibnds, 1, 2) + 1); 1];
u = u / prod(n);

% Now we calculate the A matrix.
% Let mid(i) = (hi(i) + lo(i))/2.  Let midv(i) = v_i(lo(i)) + v_i(hi(i))/2.
% The off-diagonal elements of A are easy:
%     SUM_{ii=lo(i)}^hi(i) SUM_{jj=lo(j)}^hi(j) v_i(ii)*v_j(jj)
%   = (SUM_{ii=lo(i)}^hi(i) v_i(ii)) * (SUM_{jj=lo(j)}^hi(j) v_j(jj))
%   = n(i)*midv(i) * n(j)*midv(j).
% We multiply this by prod_{k NotIn {i,j}} n(k), to get the total number
% of such elements, and then divide by prod(n) to normalize.  Thus, the
% off-diagonal elements of A are given by:
%     midv * midv'.
midv = sum(vbnds,2)/2;
A = midv * midv';

% The diagonal elements of A are given by
%     SoM_{j=lo(i)}^{hi(i)} v_i(i)^2
%   = SUM_{j=lo(i)}^{hi(i)} (v0(i) + i*dv(i))^2
%   = SUM_{j=lo(i)}^{hi(i)} (v0(i)^2 + 2*v0(i)*dv(i)*j + (dv(i)*j)^2)
%   = n(i)*v0(i)^2 + 2*v0(i)*dv(i)*(SUM_{j=lo(i)}^{hi(i)} j)
%                 + dv(i)^2 * (SUM_{j=lo(i)}^{hi(i)} j^2).
%   = n(i)*v0(i)^2 + 2*n(i)*v0(i)*dv(i)*mid(i)
%       + dv(i)^2 * (   2*(hi(i)^3 - lo(i)^3) + 3*(hi(i)^2 + lo(i)^2)
%      + (hi(i) - lo(i))
%          )/6
%  Again, the last row is a special case.  We just put a 1 there.
n   = n(1:end-1);
mid = sum(ibnds,2)/2;
nn  = size(A,1);
dd  = M.GRID.v0 .* (M.GRID.v0 + (2 .* M.GRID.dv .* mid)) ...
    + M.GRID.dv.^2 .* (diff(2*ibnds.^3 + ibnds, 1, 2) + 3*sum(ibnds.^2, 2)) ...
    ./(6*n);
A(1:(nn+1):nn^2) = [dd', 1];

% Finally, calculate b
% singleton dimensions cause singularities in A, remove them
d1 = ibnds(:,1) == ibnds(:,2);  % the singular dimensions.
d2 = [find(ibnds(:,1) ~= ibnds(:,2)); nv+1]; % the non-singular ones
b = zeros(size(u));
b(d2) = A(d2,d2) \ u(d2);
b(d1) = 0;

% end of linfit


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lft_simu_fit_ycv_help(YC, ibnds, i, j, nv, k)
%  returns a vector which is
%    SUM_{i1=lo(1)}^hi(1) SUM_{i2=lo(2)}^hi(2)... SUM_{im=lo(m)}^hi(m) (
%      y(i1, i2, .. im)
%  The sums skip dimension j, and instead return a vector corresponding
%  to values elements lo(j):hi(j) of the sum over the other dimensions.
%  We calculate the sum by using the cumulative values of YC at the
%  corners of the hypercube described by ibnds.  This involves adding
%  or subtracting the value a each corner to the final sum.  The corner
%  that corresponds to maximizing all variables (the top corner) is added.
%  Other corners are added if a path from the top corner to the current
%  corner has even length, and subtracted if the length of the path is odd.
function ycv = lft_simu_fit_ycv_help(YC, ibnds, i, j, nv, k)
if(k == j)
  ycv = lft_simu_fit_ycv_help(YC, ibnds, i, j, nv, k+1);
  % by chaoyan, i(k)=0 when ibnds(k,1)=0
  if(ibnds(k,1) == 0) 
    ycv = [0; ycv]; 
    assert(false); % should never be called.
  end;
elseif(k <= nv)
  i{k} = ibnds(k,1)-1;
  if(i{k} == 0)
    ycv1 = zeros(size(i{j}))';
  else
    ycv1 = lft_simu_fit_ycv_help(YC, ibnds, i, j, nv, k+1);
  end;
  i{k} = ibnds(k,2);
  ycv2 = lft_simu_fit_ycv_help(YC, ibnds, i, j, nv, k+1);
  ycv = ycv2 - ycv1;
else
  ycv = reshape(YC(i{:}), [], 1);
end;
% end of lft_simu_fit_ycv_help
