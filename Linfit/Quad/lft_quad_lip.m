function [b,err] = lft_quad_lip(model,bbox,lp,useErr,tol)
% [b,err] = lft_quad_lip(model,bbox,lp,useErr,tol)
% The function computes a linear approximation with minimized L2 norm error.
% The N-n quadratic function to be approximated is:
%  f(x1,...,xn) = sum(sum( model.data(rnd(x1),...,rnd(xn)) .*...
%              (dx1,...dxn,1)*(dx1;...;dxn,1) ))
%  where rnd(xi) = round((xi-model.v0(i))/model.dv(i))
%  and dxi = (xi-model.v0(i))/model.dv(i) - rnd(xi), -0.5 <= dxi < 0.5
%
% That is, for a point (x1,...xn), it is normalized to the modle grid firstly
% then find the nearest grid point and use its corresponding quadratic parameters.
%
% B = model.data is a N-n cells array, each component is a (n+1)x(n+1) matrix
% However, we use only the lower triangular because it is symmetric matrix.(sum bij and bji)
% For example, the quadratic function in 3D is
%  [b1,...,b10] * [x1^2; x1x2; x1x3; x1; x2^2; x2x3; x2; x3^2; x3; 1]
%
% The linear approxmation is:
%  ff = b'[x1,...,xn,1] \pm err
%
% Input
%  model:  N-n quadratic function (usually from ccm_getModel function). requires fields
%            model.data stores the quadratic parameters of each grid
%            model.err stores error terms used when useErr is true. 
%            model.GRID.v0 is the start value of each dimension (nx1 vector)
%            model.GRID.dv is the grid size of each diemension (nx1 vector)
%            model.GRID.nv is the number of grids of each dimension (nx1 vector)
%  bbox:  nx2 matrix. bbox(i,1) is the lower bound of ith dimension and bbox(i,2) is the upper bound.
%         bbox defines a hyper-rectangle, within which we compute the approximation
%  lp:    constraint of region
%  useErr:add the quadratic arpproximation
%  tol:   The maximum gap of upper and lower bounds of error function. 5\% difference by default
%
% Output
%  b:    approximate the quadratic function within bbox as b'*[x1,...xn,1]
%  err:  The maximum error (f-ff or ff-f, because it is balanced)
%

% check parameters
if(nargin<2)
    error('not enough parameters');
end
if(any(bbox(:,2)<bbox(:,1)))
    error('infeasible region');
end
if(nargin<3)
    lp = [];
end
if(nargin<4||isempty(useErr))
  useErr = false;
end
if(nargin<5||isempty(tol))
    tol = 0.05;
end
% Integrator is zero when xl = xh. We bloat it a little to solve the problem
eps = 1e-3;
ind = find(bbox(:,2)<=bbox(:,1)+eps);
bbox(ind,:) = bbox(ind,:)+repmat([-eps,eps],length(ind),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute dimension or model related data
n = length(model.GRID.nv);
k = (n+1)*(n+2)/2; % # of quadratic parameters
qpind = tril(true(n+1,n+1)); % lower triangular elements, the order save in model.data.

% find all cubes covered by bbox
grid = bbox ./ [model.GRID.dv,model.GRID.dv] - round([model.GRID.v0,model.GRID.v0]./[model.GRID.dv,model.GRID.dv]); %round off
bboxm = [round(grid(:,1)), round(grid(:,2))-(mod(grid(:,2),1)==0.5)]; % avoid empty cube because round(0.5)=1
bboxm = max(1,min(repmat(model.GRID.nv,1,2),bboxm));
nz = (bboxm(:,2)-bboxm(:,1))+1;
nc = prod(nz); % # of cubes

% computes lower bound, upper bound and center of each cube.
L = zeros(n,nc);
H = zeros(n,nc);
S = zeros(n,nc); % center of cube
ind = cell(n,1);
for i=1:n
    % integral range for each cell. ncxn matrix
    SIZ = reshape(nz,1,[]); SIZ(i) = 1;
    SIZ2 = ones(size(SIZ)); SIZ2(i) = nz(i);
    L(i,:) = reshape(repmat(reshape([grid(i,1)-bboxm(i,1),repmat(-0.5,1,nz(i)-1)],SIZ2),SIZ),1,[]);
    H(i,:) = reshape(repmat(reshape([repmat(0.5,1,nz(i)-1),grid(i,end)-bboxm(i,2)],SIZ2),SIZ),1,[]);
    S(i,:) = reshape(repmat(reshape((bboxm(i,1):bboxm(i,2))+round(model.GRID.v0(i)/model.GRID.dv(i)),SIZ2),SIZ),1,[]); % v0/dv==-1
    ind(i) = {bboxm(i,1):bboxm(i,2)};
end
% extract quadratic parameters for each cube.
Bs = reshape(model.data(ind{:}),[],1);
Bs = reshape(cell2mat(Bs),[],nc); % (kxnc) one column for a cube.
Errs = reshape(model.err(ind{:}),1,[]); %error of cube

% compute vertices for all cubes.
choice = 2; np = choice^n;
vind = mod(floor( repmat(0:(np-1),n,1)./repmat(choice.^(((n-1):-1:0)'),1,np) ), choice);
vind = sub2ind([n,choice],repmat((1:n)',1,np),vind+1);
cx = mat2cell([L,H],ones(1,n),nc*ones(1,choice));
x = cell2mat(cx(vind));

% remove cubes voilates the lp
% intersect(lp,cube)==null <==> at least one constraint is not satisfied by all vertices
% NOTE: each point is tested 2^n times (simple code). It can be optimized
if(~isempty(lp))
    ox = (x+repmat(S,1,np)).*repmat(model.GRID.dv,1,nc*np); % restore the absolute value
    ind = lp.A*ox > repmat(lp.b,1,nc*np)+1e-9; % 1e-9: voilate clearly 
    ind = all(reshape(ind,[],np),2);  % all cubes voilate a constraint
    ind = ~any(reshape(ind,[],nc),1); % any constraint is voilated
    if(any(ind))
        L = L(:,ind); H = H(:,ind);
        S = S(:,ind); Bs = Bs(:,ind); Errs = Errs(:,ind);
        x = x(:,repmat(ind,1,np));
        nc = size(L,2);
    else % intersection lp and bbox might be empty because of round-off error.
        %fprintf('intersection of lp and bbox is empty\n');
        warning('intersection of lp and bbox is empty\n');
        lp = [];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the best linear approximation to minize L2 norm error
% min e2f = \int e^2 dX where, e = f-ff
% de2f/db = 2*\int (f-ff)x_i dX = 0 ==> \int fx_i dX = \int ffx_i dX
% Therfore, we have a linear system
%  \int ffX dX = M*b
%  \int fX dX = c
%  b = M\c

% Compute M matrix
%   For the polynomial integral, we have the properties that
%    \int x1^k1...xn^kn  = \prod_{i=1}^{n} (dx_i^{ki+1}/(ki+1))
%  where df = f(hi) - f(lo)
%  We use <i1,...,in> to represent prod(dx1^i1,...,dxn^in)/prod(i1,...,in)
%   We have ff = [x1,...,xn,1]*b, let <ei> = [0,...,1,...,0] represent xi <e0> for 1 and <e1> for [1,...,1]
%    ff = [e1,e2,...,en,e0]*b
%  therefore:
%  |int ff*e1*b dX|     |e1|
%  |int ff*e2*b dX|     |e2|
%  |    ...      | =(\int |..|*[e1,e2,...,en,1] dX) * b;
%  |int ff*en*b dX|    |en|
%  |int ff*e0*b dX|     |e0|
%           =[(E*E')+1]*b;

% First, we generate this pattern
En   = [eye(n);zeros(1,n)];
EEn  = repmat(En,1,n+1);
EEn2 = repmat(reshape(En',1,[]),n+1,1);
pM   = EEn+EEn2; % add 1 for integral
cpM  = mat2cell(pM,ones(1,n+1),n*ones(1,n+1)); % (n+1)x(n+1) cell, each cell is 1xn vector.

% Second, compute the integral
if(isempty(lp)) % one integral
    diffs = [diff(bbox,[],2), diff(bbox.^2,[],2)/2, diff(bbox.^3,[],2)/3];
    dind = sub2ind([n,3],repmat(1:n,n+1,n+1),pM+1); % compute the indices of diffs used in M
    M = reshape(prod(reshape(diffs(dind)',n,[]),1),n+1,n+1)'; % prod from x_1 to x_n
else % sum over all cubes
    oL = (L+S).*repmat(model.GRID.dv,1,nc);  % cube coordinate -> real
    oH = (H+S).*repmat(model.GRID.dv,1,nc);
    diffs = [oH-oL; (oH.^2-oL.^2)/2; (oH.^3-oL.^3)/3];
    dind = sub2ind([n,3],repmat(1:n,n+1,n+1),pM+1); % indices for each cube
    dind = reshape(dind',[],1); % place integral <> together
    dind = sub2ind([3*n,nc],repmat(dind,1,nc),repmat(1:nc,(n+1)^2*n,1));
    M = diffs(dind);
    M = prod(reshape(M,n,[]),1); % prod integral
    M = sum(reshape(M,[],nc),2); % sum over all cubes
    M = reshape(M,n+1,n+1);
end

% Compute the c vector
%  Because the quadratic polynoimal is different for different cells, we have to compute
%  the integral for each cell and sum the result. For a cell, the integral is
%  In+1 = \int f dX
%    = \int sum(sum( (X*X').*B ))  dX
%     Notice: here x is the corresponding varible in the grid, dx in fact.
%        when we change variable in integral, both integral range and function should be changed
%        we should compute df/dx = dv, therefore, we move it to the last step
%    = sum(sum( (\int X*X' dX).*B ))
%    = sum(sum( (\int E*E' dX).*B ))
%    = sum(sum( [(E*E')+1].*B ))
%  However, model only store lower triangle of B, therefore, we should pickup the lower triangle only
%    Notice:  this is the same pattern used to compute matrix M
%    = sum(sum( [(E*E')+1](ind).*B ))
%   Ii = \int f Xi dX = \int f (xi+S(xi))*dv dX
%    = dv(i) * (\int f*xi dX + S(xi)*In+1);
%     |       |f*e1|        |f*e0|    |    |dv(1)|
%     |     |f*e2|        |f*e0|    |    |dv(2)|
%  I =| (\int |... | dX) + (\int |... | dX).*S | .* | ... |
%     |     |f*en|        |f*e0|    |    |dv(n)|
%     |     |f*e0|        |f*e0|    |    |dv(0)| = 1
%    =sum(sum( [(E*E'+1).*E +(E*E'+1).*S].*B))*dv

% First, generate the pattern for all cubes
pc0 = cpM(qpind); % pickup the lower triangle elements. kx1 cell;
pc0 = cell2mat(pc0'); % pattern for c_{n+1}
pc = [repmat(pc0,n,1)+repmat(eye(n),1,k); pc0]; % pattern for c, add e_i.

% Second, compute the integral
diffs = [H-L; (H.^2-L.^2)/2; (H.^3-L.^3)/3; (H.^4-L.^4)/4]; % [dX; dX^2/4; dX^3/3; dX^4/4]
dind = sub2ind([n,4],repmat(1:n,n+1,k),pc+1); % indices of within each column (cube)
dind = reshape(dind',[],1); % tranpose to place integral <> together
dind = sub2ind([4*n,nc],repmat(dind,1,nc),repmat(1:nc,(n+1)*k*n,1));% indices of all cubes, column for cube

I = diffs(dind); % pickup all integral terms
I = prod(reshape(I,n,[]),1); % compute the prod from x_1 to x_n
I = reshape(I,[],nc); % column for one cube

% Third, compute c
c = repmat(Bs,n+1,1).*I; % times by parameters
c = reshape(sum(reshape(c,k,[]),1),[],nc); % sum over all k polynomial terms
c(1:n,:) = c(1:n,:) + repmat(c(end,:),n,1).*S;  % add c_{n+1} to c_i
c = sum(c,2); % sum over all cubes

c = c.*([model.GRID.dv;1]); % cube length
c = c.*prod(model.GRID.dv); % constant term because of changing variable

% Finally, find the best fit
b = M\c;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the error
% f(x) = u'*B*u
% ff(x) = b'[xx;1]; xx = (x+s).*model.GRID.dv
% e = f-ff = u'*eA*u;

% map indices from B to full version \hat{B};
map = zeros(n+1,n+1);
map(qpind) = 1:k;  map = map+map'-diag(diag(map));
cind = map(:,end);

% symmetric matrix eA
eBs = Bs;
xdv = repmat(b(1:end-1).*model.GRID.dv,1,nc); % b_i*dv_i
eBs(cind(1:end-1),:) = eBs(cind(1:end-1),:) - xdv; % -b_i*x_i*dv_i
eBs(cind(end),:) = eBs(cind(end),:) - (b(end) + sum(xdv.*S,1)); % constant term
eAs = repmat(reshape((1+eye(n+1))/2,[],1),1,nc).*eBs(map,:); % make it symmetric

% evaluate error for all grid points
% Notice: e is not continous over cubes because cosine window is ignored.
u = [x; ones(1,nc*np)];
dv = repmat(u,n+1,1).*reshape(repmat(reshape(u,[],1),1,n+1)',(n+1)^2,[]);
lerr = sum(repmat(eAs,1,np).*dv,1);
lerr = [min(lerr);max(lerr)];

counter = 1;
while(true)
    x0 = (L+H)/2; ds = (H-L)/2;

    % compute lower bound
    u = [x0;ones(1,nc)];
    dv = repmat(u,n+1,1).*reshape(repmat(reshape(u,[],1),1,n+1)',(n+1)^2,[]);
    e0 = sum(eAs.*dv,1);
    lerr = [ min(lerr(1),min(e0,[],2)); max(lerr(2),max(e0,[],2)) ];

    % compute upper bound
    [uemin,uemax] = quadBounds(eAs,x0,ds,e0);
    uerr = [min(uemin);max(uemax)];

    r = diff(uerr)/diff(lerr);
    if(r>1+tol)
        if(counter>5 || nc >1e6)
            [b,err] = lft_quad_pt(model,bbox,lp);
            return;
        end
        % remove cube with smaller error than lower bound
        ind = (uemin <lerr(1) | uemax >lerr(2));
        L = L(:,ind); %H = H(:,ind);
        x0 = x0(:,ind); ds = ds(:,ind);
        eAs = eAs(:,ind); S = S(:,ind);
        nc = size(L,2);

        % split each cube to 2^n small ones
        newL = mat2cell([L,x0],ones(1,n),nc*ones(1,choice));
        L = cell2mat(newL(vind));
        H = L+repmat(ds,1,np);
        S = repmat(S,1,np);
        eAs = repmat(eAs,1,np);
        nc = nc*np;

        % remove cube voilates lp
        if(false && ~isempty(lp))
            cx = mat2cell([L,H],ones(1,n),nc*ones(1,choice));
            x = cell2mat(cx(vind));
            ox = (x+repmat(S,1,np)).*repmat(model.GRID.dv,1,nc*np); % restore the absolute value
            ind = lp.A*ox > repmat(lp.b,1,nc*np);
            ind = all(reshape(ind,[],np),2); % all cubes voilates a constraint
            ind = ~any(reshape(ind,[],nc),1); % any constraint is not satisfied by all vertices
            if(isempty(find(ind,1))) % all cube with larger error than lerr does not satisfies the lp
                err = lerr;
                break;
            end
            L = L(:,ind); H = H(:,ind);
            S = S(:,ind); eAs = eAs(:,ind);
            nc = size(L,2);
        end
    else
        err = uerr; % over approximated
        break;
    end
    counter = counter+1;
end

% Make the error balance.
b(end) = b(end)+(mean(err));
err = diff(err)/2;
if(useErr)
  err = err+max(Errs);
end
