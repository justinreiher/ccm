function [c,err] = lft_brock_lls(brockett,ibnds,region)
% [c,err] = lft_brock_lls(brockett,ibnds,region)
% Find opitmal L2 norm error by least square method based on polygon representation.
switch(region)
  case {1,3}
    upl = brockett.get('oupl');
    lpl = brockett.get('olpl');
  case 2 
    upl = brockett.get('oupl'); 
    lpl = brockett.get('iupl');
  case 4
    upl = brockett.get('olpl'); 
    lpl = brockett.get('ilpl');
  otherwise 
    error('region must be 1-4');
end

% crop polyline
upl = polyline_crop(upl,ibnds);
lpl = polyline_crop(lpl,ibnds);

% find all vertices
x = unique([upl(1,:),lpl(1,:)]); % all x points
% NOTE: x is scalar when ibnds(1) = ibnds(2)
if(numel(x)==1) 
  yu = upl(2,1);
  yl = lpl(2,1);
else 
  yu = interp1(upl(1,:),upl(2,:),x); 
  yl = interp1(lpl(1,:),lpl(2,:),x);
end

% find linear coefficient using the middle value
ym = (yu+yl)/2;
c = polyfit(x,ym,1);
% NOTE c can be Inf if diff(x) is around eps
if(any(isnan(c))||any(isinf(c)))
  c = [0,0];
end

% compute error
yy = [yu,yl]; xx = [x,x];
diffs = yy - c*[xx;ones(1,length(xx))];
err = [min(diffs),max(diffs)];

% shift constant term 
c(2) = c(2)+mean(err);
err = diff(err)/2;
