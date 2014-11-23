function test_lft_quad
  disp('This function test lft_quad functions');

  names = {'nmos','pmos','inv_2','snmos','spmos'};
  ccm_cfg('set','type','quad'); 
  for i=1:length(names)
    models{i} = ccm_getModel(names{i});
  end  

  disp('Check lft_quad function to make sure it includes the ids value');
  test_lft_quad_ids(models);
  disp('Compare the accuracy and performance of pt and lip methods. This may take a while'); 
  test_lft_quad_method(models);
  disp('Compare the error of with/without quadratic modeling error. This may take a while'); 
  test_lft_quad_err(models);

function test_lft_quad_ids(models)
  % this function compare inclusion and ids function
  methods = {'pt','lip'};
  useErr = false;
  for mid = 1:length(methods)
    method = methods{mid};
    disp([' checking method ',method]);
    for i=1:length(models) 
      m = models{i}; nn = length(m.GRID.v0); 
      for j=1:nn
        pt = test_lft_grid_randPt(m.GRID); b  = test_lft_grid_randBox(m.GRID);
        bbox = repmat(pt,1,2); bbox(j,:) = b(j,:); 
        [c,err] = lft_quad(m,bbox,method,[],useErr);
        grids = b(j,1):1e-4:b(j,2); ng=length(grids);
        V = repmat(pt,1,ng); V(j,:) = grids;
        ids1 = interp_data(m,V,'lookup')'; % same model with lft_quad
        ids2 = interp_data(m,V,'coswin')';
  
        % check result
        cpt = c'*[V;ones(1,ng)]; upt = cpt+err; lpt = cpt-err;
        if(any(ids1>upt)||any(ids1<lpt))
          error('The inclusion does not contain the ids function');
        end
        % plot result
        figure; hold on; grid on; 
        plot(grids,ids1,'b'); plot(grids,ids2,'c');
        plot(grids,cpt,'g'); plot(grids,lpt,'r'); plot(grids,upt,'k');
        xlabel(['v',num2str(j)]); ylabel('ids');
        legend('ids-lookup','ids-coswin','inclusion-center','inclusion-lower','inclusion-upper');
      end
    end
  end

function test_lft_quad_method(models)
  % this function compare lip and pt methods
  methods = {'pt','lip'}; N = 1000;
  fig1 = figure; hold on; grid on; xlabel('calls'); ylabel('ratio of error');  
  fig2 = figure; hold on; grid on; xlabel('models'); ylabel('ratio of time');
  for i=1:length(models) 
    m = models{i}; 
    bboxes = cell(N,1); err1 = zeros(N,1); err2 = err1;
    for j=1:N
      bboxes{j} = test_lft_grid_randBox(m.GRID); 
    end 
    t = cputime;
    for j=1:N 
      [c,err1(j,1)] = lft_quad(m,bboxes{j},methods{1}); 
    end
    t1 = cputime - t;
    t = cputime;
    for j=1:N 
      [c,err2(j,1)] = lft_quad(m,bboxes{j},methods{2}); 
    end
    t2 = cputime -t;
    figure(fig1); plot(1:N,err2./err1,utils_palette(i));
    figure(fig2); plot(i,t2/t1,'r*'); 
  end


function test_lft_quad_err(models)
  % this function compare function with or without err.
  methods = {'pt','lip'}; N = 1000;
  fig1 = figure; hold on; grid on; xlabel('calls'); ylabel('ratio of error');  title(methods{1});
  fig2 = figure; hold on; grid on; xlabel('calls'); ylabel('ratio of error');  title(methods{2});
  for i=1:length(models) 
    m = models{i}; 
    bboxes = cell(N,1); err1 = zeros(N,1); err2 = err1;
    for j=1:N
      bboxes{j} = test_lft_grid_randBox(m.GRID); 
    end 
    for j=1:N 
      [c,err1(j,1)] = lft_quad(m,bboxes{j},methods{1}); 
      [c,err2(j,1)] = lft_quad(m,bboxes{j},methods{1},[],true); 
    end
    figure(fig1); plot(1:N,sort(log(err2./err1)),utils_palette(i)); 
    for j=1:N 
      [c,err1(j,1)] = lft_quad(m,bboxes{j},methods{2}); 
      [c,err2(j,1)] = lft_quad(m,bboxes{j},methods{2},[],true); 
    end
    figure(fig2); plot(1:N,sort(log(err2./err1)),utils_palette(i)); 
  end

