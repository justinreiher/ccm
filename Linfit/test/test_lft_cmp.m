function test_lft_cmp
  disp('This function compares lft_simu and lft_quad functions');

  names = {'nmos','pmos','inv_2','snmos','spmos'};
  ccm_cfg('set','type','simu'); 
  for i=1:length(names)
    models{1,i} = ccm_getModel(names{i});
  end  
  ccm_cfg('set','type','quad'); 
  for i=1:length(names)
    models{2,i} = ccm_getModel(names{i});
  end  

  test_lft_cmp_method(models);


function test_lft_cmp_method(models)
  N = 1000;
  fig1 = figure; hold on; grid on; xlabel('calls'); ylabel('ratio of error');  
  fig2 = figure; hold on; grid on; xlabel('models'); ylabel('ratio of time');
  for i=1:size(models,2) 
    m1 = models{1,i}; m2 = models{2,i}; 
    bboxes = cell(N,1); err1 = zeros(N,1); err2 = err1;
    for j=1:N
      bboxes{j} = test_lft_grid_randBox(m1.GRID); 
    end 
    t = cputime;
    for j=1:N 
      [c,err1(j,1)] = lft_simu(m1,bboxes{j});
    end
    t1 = cputime - t;
    t = cputime;
    for j=1:N 
      [c,err2(j,1)] = lft_quad(m2,bboxes{j});
    end
    t2 = cputime -t;
    figure(fig1); plot(1:N,log(sort(err2./err1)),utils_palette(i));
    figure(fig2); plot(i,t2/t1,'r*'); 
  end

