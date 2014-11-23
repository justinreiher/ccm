function lib_coho_test
  % this may take 10 mins
  assert(strcmp(ccm_cfg('get','process'),'180nm'));
  disp('=====================================================================');
  test_source
  disp('=====================================================================');
  test_leaf
  disp('=====================================================================');
  test_component
  disp('=====================================================================');
  test_circuit
end

% voltage source and current source 
function test_source
  vnn = 1.8;

  disp('Test voltage source');
  c = inverter('inv',1e-5);

  disp('  test vsrc');
  p = vsrc('input',vnn); 
  s = testbench(c,{c.i},{p});
  [t,v,dv] = s.simulate([0,1e-8]);
  figure; hold on; 
  plot(t,v(:,1)); plot(t,v(:,2),'r');

  disp('  test vpulse');
  p = vpulse('input',[0,vnn],[2,8,2,8]*1e-10,1e-10); 
  s = testbench(c,{c.i},{p});
  [t,v,dv] = s.simulate([0,1e-8]);
  figure; hold on;
  plot(t,v(:,1)); plot(t,v(:,2),'r');

  disp('  test vsin');
  i = vsin('input',1e10,pi/2,vnn/2,1); 
  s = testbench(c,{c.i},{i});
  [t,v,dv] = s.simulate([0,1e-8]);
  figure; hold on;
  plot(t,v(:,1)); plot(t,v(:,2),'r');

  disp('  test Brockett Annulus');
  a = annulus('input',[1e10,2e10],[0,0.15,vnn-0.15,vnn],1e-10);
  a.rand(1e-8);
  s = testbench(c,{c.i},{a});
  [t,v,dv] = s.simulate([0,1e-8]);
  figure; hold on;
  plot(t,v(:,1)); plot(t,v(:,2),'r');
  figure; hold on;
  a.display;
  plot(v(:,1),dv(:,1),'r');
  r = coho_shape([0,0.2;0,0.2],[]);
  [c,err] = s.dV_ldi(r,1);

  disp('  voltage sources could be two terminal, but testbench does not support it currently')
  

  disp('Test current source');
  isrc = lib_coho_test_circuit_isrc('lib_coho_test_circuit_isrc',1e-3,1e-9); 
  isrct = testbench(isrc,{},{});
  [t,v] = isrct.simulate([0,1e-5]);
  figure; plot(t,v);

  disp('  current sources could be two terminals');
  isrc2 = lib_coho_test_circuit_isrc2('lib_coho_test_circuit_isrc2',1e-3,1e-9); 
  isrct2 = testbench(isrc2,{},{});
  [t,v] = isrct2.simulate([0,1e-5]);
  figure; plot(t,v);

end

function test_leaf
  vnn = 1.8;
  vdd = vsrc('vdd',vnn); gnd = vsrc('gnd',0);
  i = vpulse('input',[0,vnn],[2,8,2,8]*1e-10,1e-10); 

  % tsmc(simu), tsmc(quad), ptm(simu), ptm(quad)
  %    10         1e-8        10        1e-8
  %   osc off 0,  osc on 0  osc off 0   negative (wrong?
  % because ptm(quad) model give i=0(slightly positive) even vd is negative.
  % This is because the range for quad is [0,1.8], fixed by widdening to [-0.1,1.9].
  disp('Test nmos');
  c = nmos('n',1e-5);
  s = testbench(c,{c.s,c.g},{gnd,vdd});
  [t,v] = s.simulate([0,1e-10]);
  figure;plot(t,v(:,3),'r');
  
  disp('Test pmos');
  c = pmos('p',1e-5);
  s = testbench(c,{c.s,c.g},{vdd,gnd});
  [t,v] = s.simulate([0,1e-10]);
  figure;plot(t,v(:,3),'r');

  disp('Test inverter');
  c = inverter('inv',1e-5);
  s = testbench(c,{c.i},{i}); %TODO: not cell if only one input?
  [t,v] = s.simulate([0,4e-9]);
  figure; plot(t,v);

  disp('Test capacitor');
  c = capacitor('cap',1e-12);
  s = testbench(c,{c.x},{vdd}); 
  [t,v] = s.simulate([0,1e-9]);
  assert(all(v==vnn))

  disp('Test resistor');
  r = resistor('res',1e-6);
  s = testbench(c,{c.x},{vdd}); 
  [t,v] = s.simulate([0,1e-9]);
  assert(all(v==vnn))

  disp('Test inductor. It is not supported now');

  disp('Test vecotorization');
  c = nmos('nmos',1e-5); c2 = nmos('nmos',1e-5);
  i = c.I(rand(3,100));
  c.I_objs({c,c2},rand(3,100,2));
end

function test_component
  vnn = 1.8;
  vdd = vsrc('vdd',vnn); gnd = vsrc('gnd',0);
  % make sure the clock is not faster than the delay
  i1 = vpulse('input',[0,vnn],[1,8,1,8]*1e-10); 
  i2 = vpulse('input',[0,vnn],[0.5,4.5,0.5,4.5]*1e-10,1e-10); 

  disp('Test INV'); 
  c = INV('inv',1e-5);
  s = testbench(c,{c.i},{i1}); 
  [t,v] = s.simulate([0,1e-9]);
  c2 = inverter('inv',1e-5);
  s2 = testbench(c2,{c2.i},{i1}); 
  [t2,v2] = s2.simulate([0,2e-9]);
  % could be different: 1) sum(intep()) v.s. interp(sum()). 
  %assert(all(t==t2)&all(v(:)==v2(:))); 
  figure; hold on; plot(t,v(:,2)); plot(t2,v2(:,2),'r');

  tic
  disp('Test NAND. This may take 1 minute'); 
  c = NAND('nand',1e-5);
  s = testbench(c,{c.i1,c.i2},{i1,i2}); 
  [t,v] = s.simulate([0,2e-9]);
  figure; plot(t,v);
  toc

  tic
  disp('Test NOR. This may take 1 minute'); 
  c = NOR('nor',1e-5);
  s = testbench(c,{c.i1,c.i2},{i1,i2}); 
  [t,v] = s.simulate([0,2e-9]);
  figure; plot(t,v);
  toc
  
  tic
  disp('Test LATCH. This may take 3 minutes'); 
  c = LATCH('latch',1e-5);
  s = testbench(c,{c.d,c.en},{i1,i2}); 
  [t,v] = s.simulate([0,2e-9]);
  inds = c.find_port_index({c.d,c.en,c.q});
  figure; plot(t,v(:,inds));
  toc

  % NOTE: compared with latch, the time step is 2x smaller 
  tic
  disp('Test FF. This may take 5 minutes'); 
  c = FF('latch',1e-5);
  s = testbench(c,{c.d,c.clk},{i1,i2}); 
  [t,v] = s.simulate([0,2e-9]);
  inds = c.find_port_index({c.d,c.clk,c.q});
  figure; plot(t,v(:,inds));
  toc

  tic
  disp('Test CELEMENT. This may take 1 minutes'); 
  c = CELEMENT('c',1e-5);
  s = testbench(c,{c.i1,c.i2},{i1,i2}); 
  [t,v] = s.simulate([0,2e-9]);
  inds = c.find_port_index({c.i1,c.i2,c.o});
  figure; plot(t,v(:,inds));
  toc
end

function test_circuit
  disp('Test RRO');
  N = 20;
  c = RRO('rro',N,1e-5);
  s = testbench(c,{},{});
  v0 = rand(2*N,1)*1e-3;
  [t,v] = s.simulate([0,1e-8],v0);
  figure; hold on; plot(t,v(:,1)); plot(t,v(:,end),'r');

  disp('Test IRO');
  c = IRO('rro',5,1e-5);
  s = testbench(c,{},{});
  v0 = rand(5,1)*1e-3;
  [t,v] = s.simulate([0,2e-9],v0);
  figure; hold on; plot(t,v(:,1)); plot(t,v(:,end),'r');
end

