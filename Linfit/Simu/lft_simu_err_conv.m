function [err, pt] = lft_simu_err_conv(M, u, b)
% [err,pt] = lft_simu_err_conv(M,u,b)
%   This function compute the linearization error.
%       err(1) = min_{i IN u} M.data(i) - b'* v(i)
%       err(2) = max_{i IN u} M.data(i) - b'* v(i)
%       pt is the point that reach the extrme error.
%   This function is more efficient than lft_simu_err_bf when u is large.
%   However, lft_simu_err_bf is more efficient when u is small.
%   It works only for MOS ids with assumption of convexity of vg vd and vs.
[err, pt] = lft_err_conv_help(M,u,b);
err = -[err(2);err(1)]; pt=[pt(2,:);pt(1,:)];

function [err, pt] = lft_err_conv_help(M, u, b)
% err(1) = min_{i IN u} b'*v(i) - M.data(i)
% err(2) = max_{i IN u} b'*v(i) - M.data(i)
%
% BUG: chaoyan, the method might be sensitive to M.GRID.dv or other round off error.
% The error is small now. Try s=0,g=0-1.8,d=1.8. the approximation is not correct.(lft_simu_err_bf) 
%
% I've run out of time.  I won't worry about generality for now.
% Instead, I'll make the following assumptions:
%  * The model has three free variables, vs, vg, and vd (in that order).
%  * ids is convex (i.e. concave upward) in vs.
%  * ids is concave (ie. concave downward) in vd.
%  * if vs < vd, then ids is convex in vg
%  * elseif vs > vd, then ids is concave in vg
%  * else vs == vd and ids is zero.
% These convexity assumptions means that the extremal values are either
% at the endpoints of the vs, vg, or vd intervals, or can be found by
% bisection.  We enumerate the endpoint combinations when needed, and
% use bisection on what's left.


% First, we find the indices for the cube.
% We take floor or ceil to make sure that our subcube of the sample
% contains the cube described by v (the components of v are continuous).
% ibnds(i,1) is the lowest index for the i^th dimension of M.YC for our cube.
% ibnds(i,2) is the highest index for the i^th dimension of M.YC.
% vbnds is the approximation of v corresponding to ibnds.
v0 = [M.GRID.v0, M.GRID.v0];  dv = [M.GRID.dv, M.GRID.dv];
ibnds = (u - v0) ./ dv;
ibnds = [floor(ibnds(:,1)), ceil(ibnds(:,2))];
%vbnds = [ v0 + dv .* ibnds; [1,1]];


% First, find the min of z = predicted - actual
%  z is concave in vs and convex in vd;
%  z is concave in vg for vs < vd and
%  z is convex in vg for vs > vd and
err = zeros(2,1); pt = zeros(2,3);
g = 2;
for s = [1,3]
  % When s=1, we calculate min b'*v(i) - M.data(i)
  % When s=3, we swap the source and drain, and repeat the min calculation.
  %   negating the result gives us max b'*v(i) - M.data(i)
  d = 4-s;
  minpt = []; minerr = [];
  for loops = 1:2
    if(loops == 1), is = ibnds(s,1);
    else is = ibnds(s,2);
    end;
    % consider the region where vs <= vd
    lo_d = min(max(floor(((M.GRID.v0(s)+M.GRID.dv(s)*is)-M.GRID.v0(d)) / dv(d)), ...
      ibnds(d,1)), ibnds(d,2));
    hi_d = ibnds(d,2);
    for loopg = 1:2
      if(loopg == 1), ig = ibnds(g,1);
      else ig = ibnds(g,2);
      end;
      if(lo_d <= hi_d)
        [min_e, min_i] = ...
          lft_simu_err_bisect([lo_d, hi_d], @(id)(lft_simu_err_eval(M, b, [is, ig, id], s)));
        if(isempty(minerr) || (min_e < minerr))
          minpt = [is, ig, min_i]; minerr = min_e;
        end;
      end;
    end; % for loopg

    % consider the region where vs >= vd
    % pick an initial value for ig
    ig = floor(sum(ibnds(g,:),2)/2);  % we'll use the midpoint (why not?)
    min_i = []; min_e = [];
    % lft_simu_err_bisect on the drain voltage to minize the function
    ranged = [ ibnds(d,1), ...
      max(min(ceil(((M.GRID.v0(s)+M.GRID.dv(s)*is)-M.GRID.v0(d)) / dv(d)), ...
      ibnds(d,2)), ibnds(d,1));
      ];
    if(ranged(1) < ranged(2))
      while(1)
        [e, id] = lft_simu_err_bisect(ranged, @(id)(lft_simu_err_eval(M, b, [is, ig, id], s)));
        [e, ig] = lft_simu_err_bisect(ibnds(2,:), @(ig)(lft_simu_err_eval(M, b, [is, ig, id], s)));
        if(isempty(min_i) || (e < min_e))
          min_i = [is, ig, id]; min_e = e;
        else break;
        end;
      end;
      if(e < minerr)
        minpt = min_i; minerr = e;
      end;
    end;
  end; % for loops
  err((s+1)/2) = minerr*(2-s);
  pt((s+1)/2,:) = [minpt(s), minpt(g), minpt(d)];
end;  % for s
% end linfitError

function ii = lft_simu_err_eval(M, b, p, s)
p = [p(s), p(2), p(4-s)];
i_est = b' * [M.GRID.v0 + p'.*M.GRID.dv; 1];
i_m   = M.data(p(1), p(2), p(3));
ii = (i_est - i_m)*(2-s);
% end ids_eval;



% [min_v, min_i] = lft_simu_err_bisect(range, f):
%  Search a convex function for its minimum value.
%  range(1) is the lower bound (inclusive) of the range to search.
%  range(2) is the upper bound (inclusive) of the range to search.
%  f is the function.
%  min_i is the value of the argument to f that minimizes f(i).
%  min_v = f(min_i).
%  The argument to f can only take on integer values.
%  range(1) and range(2) must be integers.
function [min_v, min_i] = lft_simu_err_bisect(range, f)
lo = range(1); hi = range(2);
if(lo > hi)
  min_i = []; min_v = [];
  return;
end;
vlo = f(lo); vhi = f(hi);
if(vlo < vhi), min_i = lo; min_v = vlo;
else min_i = hi; min_v = vhi;
end;
mid = [];
while(lo < hi-1)
  if(isempty(mid))
    mid = floor((lo + hi)/2);
    vmid = f(mid);
  end;
  if((vlo <= vmid) && (vmid <= vhi))
    hi = mid; vhi = vmid; mid = [];
  elseif((vlo >= vmid) && (vmid >= vhi))
    lo = mid; vlo = vmid; mid = [];
  elseif((vmid <= vlo) && (vmid <= vhi))
    min_i = mid; min_v = vmid;
    if(mid == hi-1), return; end;
    m34 = floor((mid + hi)/2);
    v34 = f(m34);
    if(v34 < vmid),  lo = mid;   vlo = vmid; 
    else hi = m34;   vhi = v34;
    end;
    mid = [];
  else % wrong convexity
    % This might happen because the SPICE data has small, local,
    % non-convexities.  We'll test all points from lo to hi.
    % Later, I should add a test to make sure there are not too
    % many such points.
    for i = (lo+1):(hi-1)
      if(i == mid), continue; end;
      v = f(i);
      if(v < min_v), min_i = i; min_v = v; end;
    end; % for i
    return;
  end; % else, wrong convexity
end; % while
% end lft_simu_err_bisect
