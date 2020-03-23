function test_interp
  disp('This function is to test the interp_data function or circuit models'); 
  disp('Please check the result by plotted figures');

  % Load all models
  disp('loading models ...');
  names = {'nmos','pmos','snmos','spmos','inv_2'}; 
  ccm_cfg('set','type','simu');
  models = cell(2,length(names));
  for i=1:length(names)
    models{1,i} = ccm_getModel(names{i});
  end
  ccm_cfg('set','type','quad');
  for i=1:length(names)
    models{2,i} = ccm_getModel(names{i});
  end

  disp('test plot function');
  test_plot(models(:,1));

  disp('sweep 1d, compare different interp methods'); 
  test_1d_method(models);

  disp('sweep 1d, compare simu and quad models'); 
  test_1d_model(models);

  disp('sweep 2d'); 
  test_2d(models);
end

function test_plot(models)
% test plot function
  ids_sweep_1d(models{2,1},[0;1.8;1.8],'coswin',2);
  ids_sweep_2d(models{1,1},[0;1.8;1.8],'linear',[1,2]);
  ids_sweep_1d(models{2,1},[0;0;0],'coswin',2);
  ids_sweep_2d(models{1,1},[0;0;0],'linear',[1,2]);
end

% test interp_data function
function test_1d_method(models)
% sweep 1d (different methods)
  methods = {'coswin','linear','lookup'}; 
  for m = 1:length(models(:))
    model = models{m}; 
    gnd = model.META.gnd; vdd = model.META.vdd; 
    switch(length(model.GRID.v0)) 
      case 2 
        V = [gnd;vdd]; 
      case 3 
        V = [gnd;vdd;vdd]; 
      otherwise 
        V = rand(length(model.GRID.v0),1)*(vdd-gnd)+gnd; 
    end
    grids = gnd:0.005:vdd;
    for t = 1:length(model.GRID.v0); 
      fig=figure; hold on; grid on;
      s = sprintf('model=%s-%s,terminal=%d',model.META.type,model.META.name,t);
      title(s);
      xlabel(['v',num2str(t)]); ylabel('ids');
      colors = {'r','g','b'};
  
      for i =1:length(methods) 
        method = methods{i};
        [ids,grids] = ids_sweep_1d(model,V,method,t,grids); 
        plot(grids,ids,colors{i}); 
      end
      legend(methods{:});
      print(fig,'-deps2',['1dmethod_',s]); 
    end
  end
end % function

function test_1d_model(models)
% sweep 1d (different model)
  methods = {'coswin','linear','lookup'}; 

  for m=1:size(models,2)
    model1 = models{1,m}; model2 = models{2,m}; % simu and quad models
    gnd = model1.META.gnd; vdd = model2.META.vdd; 
    V = rand(length(model1.GRID.v0),1)*(vdd-gnd)+gnd;
    grids = gnd:0.005:vdd; 

    for t = 1:length(model1.GRID.v0);
      for i =1:length(methods) 
        method = methods{i};
        fig=figure; hold on; grid on; 
        s = sprintf('model=%s,terminal=%d,method=%s',model1.META.name,t,method);
        title(s); xlabel(['v',num2str(t)]); ylabel('ids');
        [ids1,grids1] = ids_sweep_1d(model1,V,method,t,grids); 
        [ids2,grids2] = ids_sweep_1d(model2,V,method,t,grids); 
        plot(grids1,ids1,'r'); 
        plot(grids2,ids2,'b'); 
        legend('simu','quad');
        print(fig,'-deps2',['1dmodel_',s]);
      end
    end
  end  
end % function

function test_2d(models)
% sweep 2d
  methods = {'coswin','linear','lookup'}; 

  for m = 1:length(models(:))
    model = models{m}; 
    gnd = model.META.gnd; vdd = model.META.vdd; 
    grids = {gnd:0.05:vdd,gnd:0.05:vdd};
    switch(length(model.GRID.v0)) 
      case 2
        ts = {[1,2]}; V = [gnd;vdd]; 
      case 3
        ts = {[1,2],[1,3],[2,3]}; V = [gnd;vdd;vdd]; 
      otherwise
        error('do not support now');
    end
    for tt = 1:length(ts)
      t = ts{tt};
  
      for i =1:length(methods) 
        method = methods{i};
        fig=figure; hold on; grid on;
        s = sprintf('model=%s-%s,terminal=[%d,%d],method=%s',model.META.type,model.META.name,t(1),t(2),method);
        title(s);
        xlabel(['v',num2str(t(1))]);ylabel(['v',num2str(t(2))]);zlabel('ids');
        [ids,X1,X2] = ids_sweep_2d(model,V,method,t,grids); 
        mesh(X1,X2,ids); view(-37.5,30);
        axis([min(X1(:)),max(X1(:)),min(X2(:)),max(X2(:)),min(ids(:)),max(ids(:))]);
        print(fig,'-deps2',['2d_',s]);
      end
    end
  end
end % function


function [ids,grids] = ids_sweep_1d(model,V,method,t,grids)
  if(nargin<6), grids = []; end
  func = @(model,V,method) interp_data(model,V,method); 
  if(nargout==0)
    sweep1d(model,V,method,t,grids,func);
  else
    [ids,grids] = sweep1d(model,V,method,t,grids,func);
  end
end

function [ids,X1,X2] = ids_sweep_2d(model,V,method,t,grids)
  if(nargin<6), grids = []; end
  func = @(model,V,method) interp_data(model,V,method); 
  if(nargout==0) 
    sweep2d(model,V,method,t,grids,func);
  else 
    [ids,X1,X2] = sweep2d(model,V,method,t,grids,func);
  end
end
