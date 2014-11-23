function iData = spice_device_ptm(grids,opt,tmpName)
% iData = spice_device(grids,opt,tmpName)
% This function is used to generate spice file to simulate a device
% and save its current function as a table.
  if(nargin<2||isempty(opt))
    opt = [];
  end
  if(nargin<3||isempty(tmpName))
    tmpName = 'tmp';
  end
  assert(length(grids)==3); % hspice can only sweep two terminals
  iData = zeros(length(grids{1}),length(grids{2}),length(grids{3}));
  v1 = grids{1};
  for i=1:length(v1)
    g = {v1(i),grids{2},grids{3}};
    iData(i,:,:) = spice2d(g,opt,tmpName);
    fclose('all');
  end
end % spice_device

function iData = spice2d(grids,opt,tmpName)
  str = idsStrings(grids,opt);
  iData = runSpice(str,tmpName);
end % spice2d 

function spiceData = runSpice(deckString,tmpName) 
  % open a file
  spiceName = [tmpName, '.hsp']; 
  tmpFile = fopen(spiceName, 'w'); 
  if(tmpFile < 0) 
    error(['could not open file for write: ', spiceName]); 
  end 
  % write to the file
  for i=1:length(deckString)
    fprintf(tmpFile, '%s\n', deckString{i}); 
  end
  fclose(tmpFile); 
  % run hspice 
  [status, result] = unix(['hspice ', spiceName, ' > ', tmpName, '.out']); 
  if(status ~= 0), error('hspice failed'); end 
  % load the output data 
  outName = [tmpName, '.out']; 
  spiceData = loadData(outName); 
end % runSpice 

% collect the spice data and put it into a table
function spiceData = loadData(outName) 
  spFile = fopen(outName, 'r'); 
  tline = fgets(spFile); 
  str = sscanf(tline, '%s');
  j = 1;
  while (tline ~= -1) 
    while (~strcmp(str,'x')) 
      tline = fgets(spFile); 
      if (tline == -1), return; end 
      str = sscanf(tline, '%s'); 
    end
    % skip 4 lines
    tline = fgets(spFile); tline = fgets(spFile); 
    tline = fgets(spFile); tline = fgets(spFile);
    str = sscanf(tline, '%s');
    i = 1;
    while (~strcmp(str,'y'))
      data = sscanf(tline(15:end), '%f%c');
      if (char(data(2)) == 'a') % added by chaoyan, is 'a' 1e-18? 
        scale = 1e-18;
      elseif (char(data(2)) == 'f')
        scale = 1e-15;
      elseif (char(data(2)) == 'p')
        scale = 1e-12;    
      elseif (char(data(2)) == 'n')
        scale = 1e-9;
      elseif (char(data(2)) == 'u')
        scale = 1e-6;
      elseif (char(data(2)) == 'm')
        scale = 1e-3;
      else
        scale = 1;
      end % if 
               
      % Why negative?
      spiceData(i,j) = -data(1)*scale; 
      tline = fgets(spFile); 
      str = sscanf(tline, '%s');            
      i = i+1;
    end % while
    j=j+1;
  end % while
  fclose(spFile);
end

% This function provides generate spice file to simultate a circuit 
function s = idsStrings(grids,opt)
  % process parameters
  v1 = grids{1}; v2 = grids{2}; v3 = grids{3};
  device = opt.device;
  ids = opt.ids; wid = opt.wid;

  % HSpice environment 
  preamble = { 
    ['* ',upper(device),' 018 I-V curve'],
      '',
      ['vv1 v1 0 ', num2str(v1)],
      'vv2 v2 0 0',
      'vv3 v3 0 1.8',
      '',
      '.temp 25',
      '',
      '.protect',
      '.include ''./180nm_bulk.pm'' ',
      '.unprotect',
      '' 
  };

  switch(lower(device))
    case 'nmos' 
      c = mStrings(wid, 'm1 v3 v2 v1 gnd nmos');%n mos
    case 'pmos' 
      c = mStrings(wid, 'm1 v3 v2 v1 vdd pmos');%p mos
    case 'nmos2'
      c1 = mStrings(wid(1), 'm1 vi v1 gnd gnd nmos',[1;1;0.25;0.25]);%n mos
      c2 = mStrings(wid(2), 'm2 v3 v2 vi  gnd nmos',[0.25;0.25;1;1]);%n mos
      c = [c1,c2];
    case 'pmos2'
      c1 = mStrings(wid(1), 'm1 vi v1 vdd vdd pmos',[1;1;0.25;0.25]);%p mos
      c2 = mStrings(wid(2), 'm2 v3 v2 vi  vdd pmos',[0.25;0.25;1;1]);%p mos
      c = [c1,c2];
    otherwise
      error('do not support');
  end
  
    % sweep
  if(ids) 
    vv = 'vv3'; % ids
  else 
    vv = 'vv1'; 
  end % isd
  epilog = { 
    '', 
    '.option chgtol=1e-20', 
    '.option accurate', 
    '.option dccap=1', 
    '.opt post', 
    '', 
    ['.dc vv2 ', num2str(v2(1)),' ', num2str(v2(end)),' ', num2str((v2(end)-v2(1))/(length(v2)-1)),...  
     ' vv3 ', num2str(v3(1)),' ', num2str(v3(end)), ' ',num2str((v3(end)-v3(1))/(length(v3)-1))], 
    ['.print ''ids''= I(',vv,')'], 
    '', 
    '.end', 
  }; 

  %s = {preamble{:}, c{:}, epilog{:}};
  s = [preamble; c; epilog];
end % idsStrings 

% string for capacitance
function s = mStrings(w,x,f) 
  if(nargin<3||isempty(f))
    f = ones(4,1);
  end
  fmt = '%e'; 
  first = [x, ' W =', sprintf(fmt, w), ' L=0.2u']; 
  cap = { 
    ['+ AS=', sprintf(fmt, w*0.48e-6*f(1)), ... 
      ' PS=', sprintf(fmt, (2*w+0.96e-6)*f(2))], 
    ['+ AD=', sprintf(fmt, w*0.48e-6*f(3)), ...  
      ' PD=', sprintf(fmt, (2*w+0.96e-6)*f(4))], 
    }; 
  s = [{first}; cap];
end % mStrings 
