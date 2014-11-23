function pt = test_lft_grid_randPt(GRID)
  % generate a random bbox for a GRID
  nn = length(GRID.v0);
  v0 = GRID.v0+GRID.dv;
  v1 = GRID.v0+GRID.dv.*GRID.nv;
  pt = v0+(v1-v0).*rand(nn,1);
