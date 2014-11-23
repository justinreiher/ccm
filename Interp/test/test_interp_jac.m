function test_interp_jac
  disp('This function is to test the interp_jacob* functions'); 
  disp('Please check the result by plotted figures');

  disp('loading models ...');
  names = {'nmos','pmos','snmos','spmos','inv_2'}; 
  models = cell(2,length(names));
  ccm_cfg('set','type','simu');
  for i=1:length(names)
    models{1,i} = ccm_getModel(names{i});
  end
  ccm_cfg('set','type','quad');
  for i=1:length(names)
    models{2,i} = ccm_getModel(names{i});
  end

  % compare different model
  disp('Test different models (simu vs quad)'); 
  disp('NOTE: the result for simu models are oscillating'); 
  test_model(models);

  disp('Compare interp_jac and interp_jac_num'); 
  test_method(models);

  disp('Test interp_jacob_num with all models and interp methods'); 
  test_num(models);
end

function test_model(models) 
  method = []; grids = 0:0.005:1.8;
  for i=1:size(models,2)
    m1 = models{1,i}; m2 = models{2,i};
    gnd = m1.META.gnd; vdd = m1.META.vdd;
    d = length(m1.GRID.v0);
    V = rand(d,1)*(vdd-gnd)+gnd; 
    for t = 1:d
      fig=figure; hold on; grid on; 
      s = sprintf('model=%s,terminal=%d',m1.META.name,t); 
      title(s); xlabel(['v',num2str(t)]); ylabel('jac');
      [jac1,grids1] = jac_sweep_1d(m1,V,method,t,grids);
      [jac2,grids2] = jac_sweep_1d(m2,V,method,t,grids);
      plot(grids1,jac1(t,:),'r-*'); 
      plot(grids2,jac2(t,:),'y-+'); 
      %[ids1,grids1] = ids_sweep_1d(m1,V,method,t,grids);
      %[ids2,grids2] = ids_sweep_1d(m2,V,method,t,grids);
      %plot(grids1,ids1,'b');
      %plot(grids2,ids2,'k');
      %legend('jac-simu','jac-quad','ids-simu','ids-quad');
      legend('jac-simu','jac-quad'); 
      print(fig,'-deps2',['jac1d_',s]);
    end
  end
end

function test_method(models)
  method = 'coswin'; grids = 0:0.005:1.8;
  for i=1:length(models(:))/2
    m = models{2,i}; 
    gnd = m.META.gnd; vdd = m.META.vdd;
    d = length(m.GRID.v0); 
    V = rand(d,1)*(vdd-gnd)+gnd; 
    for t = 1:d
      fig=figure; hold on; grid on; 
      s = sprintf('model=%s-%s,terminal=%d',m.META.type,m.META.name,t); 
      title(s); xlabel(['v',num2str(t)]); ylabel('jac');
      [jac1,grids1] = jac_sweep_1d(m,V,method,t,grids,false);
      [jac2,grids2] = jac_sweep_1d(m,V,method,t,grids,true);
      plot(grids1,jac1(t,:),'r-*'); 
      plot(grids2,jac2(t,:),'y-+'); 
      %[ids1,grids1] = ids_sweep_1d(m,V,method,t,grids);
      %plot(grids1,ids1,'k');
      %legend('jac-analytical','jac-numerica','ids');
      legend('jac-analytical','jac-numerica'); 
      print(fig,'-deps2',['jac1d_',s]);
    end
  end
end

% NOTE, ususally the numerical method works
%      lookup  linear coswin
% simu   x     v    x
% quad   x     v    v
% Because lookup is not continuous, linear is C^0 but not C^1.
% coswin method is C^1, but with simu model, its derivative is linear 
% combinationsin of sin function. Therefore, it always oscillates 
function test_num(models) 
  methods={'lookup','linear','coswin'}; 
  %methods={'coswin'}; % messed up to plot all
  for i=1:size(models,2)
    m1 = models{1,i}; m2 = models{2,i};
    gnd = m1.META.gnd; vdd = m1.META.vdd; 
    grids = gnd:0.005:vdd; 
    d = length(m1.GRID.v0); 
    for t=1:d 
      V = rand(d,1)*(vdd-gnd)+gnd; 
      for k=1:length(methods)
        method = methods{k};
        [jac1,grids1] = jac_sweep_1d(m1,V,method,t,grids,true);
        [jac2,grids2] = jac_sweep_1d(m2,V,method,t,grids,true);

        fig=figure; hold on; grid on; 
        s = sprintf('model=%s,method=%s,terminal=%d',m1.META.name,method,t); 
        title(s); xlabel(['v',num2str(t)]); ylabel('jac');
        plot(grids1,jac1(t,:),'r-*'); 
        plot(grids2,jac2(t,:),'y-+'); 
        legend('jac-simu','jac-quad'); 
        print(fig,'-deps2',['jac1d_',s]);
      end
    end
  end
end


function [jac,grids] = jac_sweep_1d(model,V,method,t,grids,isnum)
  if(nargin<6), grids = []; end
  if(nargin<7||isempty(isnum)), isnum = false; end
  if(isnum)
    func = @(model,V,method) interp_jacob_num(model,V,method); 
  else 
    func = @(model,V,method) interp_jacob(model,V); 
  end
  if(nargout==0)
    sweep1d(model,V,method,t,grids,func);
  else
    [jac,grids] = sweep1d(model,V,method,t,grids,func);
  end
end

% copy from test_data;
function [ids,grids] = ids_sweep_1d(model,V,method,t,grids)
  if(nargin<6), grids = []; end
  func = @(model,V,method) interp_data(model,V,method); 
  if(nargout==0)
    sweep1d(model,V,method,t,grids,func);
  else
    [ids,grids] = sweep1d(model,V,method,t,grids,func);
  end
end


