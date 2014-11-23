function bbox = test_lft_grid_randBox(GRID)
  % generate a random bbox for a GRID
  nn = length(GRID.v0);
  v0 = GRID.v0+GRID.dv;
  v1 = GRID.v0+GRID.dv.*GRID.nv;
  bbox = repmat(v0,1,2)+repmat((v1-v0),1,2).*sort(rand(nn,2),2);
