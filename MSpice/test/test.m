function test
  disp('=======================================================================');
  disp('Test testbench class'); 
  test_testbench;

  disp('=======================================================================');
  disp('Test circuit class'); 
  test_circuit;

  disp('=======================================================================');
  disp('Test flatten and vectorization');
  test_flat_vec;
  
  % NOTE: currently, circuit doesn't support source, which is only used by testbench. 
  % It's a TODO list for futhure development.
  disp('=======================================================================');
  disp('Test voltage source');
  test_src;
end

% test basic simulation functions
function test_testbench
  assert(strcmp(ccm_cfg('get','process'),'180nm'));
  vnn = 1.8; % TODO: better way to get vnn?

  disp('-----------------------------------------------------------------------');
  disp('Test simulation interface'); 

  disp('  simulate a nmos');
  c = nmos('n',1e-5); vdd = vsrc('vsrc',1.8);
  s = testbench(c,{c.g,c.d},{vdd,vdd});
  [t,v] = s.simulate([0,1e-10]);
  figure; hold on;
  plot(t,v(:,1)); plot(t,v(:,2),'r'); plot(t,v(:,3),'g');
  
  disp('  simulate an inverter'); 
  c = inverter('inv',1e-5);
  p = vpulse('input',[0,vnn],[2,8,2,8]*1e-10,1e-10); 
  s = testbench(c,{c.i},{p});
  [t,v] = s.simulate([0,1e-8]);
  figure; hold on;
  plot(t,v(:,1)); plot(t,v(:,2),'r');
  
  disp('  simulate an INV composed by a nmos and a pmos'); 
  c = INV('INV',[1e-5,2e-5]);
  p = vpulse('input',[0,vnn],[2,8,2,8]*1e-10,0); 
  s = testbench(c,{c.i},{p});
  [t,v] = s.simulate([0,1e-8]);
  figure; hold on;
  plot(t,v(:,1)); plot(t,v(:,2),'r');


  disp('-----------------------------------------------------------------------');
  disp('Test verification interface as dV_ldi');
  bbox = [0,0.2;vnn-0.2,vnn];
  r = coho_shape(bbox,[]);
  r = r.project([1,2]);

  disp('  dV_ldi of an inverter'); 
  c = inverter('inv',1e-5);
  [c,err] = c.dV_ldi(r,1);

  disp('  dV_ldi of an INV'); 
  c = INV('inv',[1e-5,2e-5]);
  a = annulus('input',[1e10,2e10],[0,0.15,vnn-0.15,vnn],1e-10);
  s = testbench(c,{c.i},{a});
  [c,err] = s.dV_ldi(r,1);
end

function test_circuit
  N = 5;
  n = nmos('n',1e-5);
  v = vsrc('vsrc',1.8); 
  i = inverter('inv',1e-5);
  iro = IRO('inv-ring-osc',N,1e-5);
  latch = LATCH('latch',1e-5);
  rro = RRO('rro',4,[1,1]*1e-5);

  disp('-----------------------------------------------------------------------');
  disp('Test interface functions');
  assert(n.ifc_simu & ~n.is_vsrc & n.ifc_verify & n.ifc_ssan)
  assert(v.ifc_simu & v.is_vsrc & ~v.ifc_verify & ~v.ifc_ssan)
  assert(rro.ifc_simu & ~rro.is_vsrc & rro.ifc_verify & rro.ifc_ssan)
  assert(latch.ifc_simu & ~latch.is_vsrc & latch.ifc_verify & latch.ifc_ssan)

  disp('-----------------------------------------------------------------------');
  disp('Test simulation interface');
  disp('  test I/C/dV functions');
  if(n.ifc_simu && ~n.is_vsrc)
    n.I(rand(3,1));
    n.C(rand(3,1));
    n.dV(rand(3,1));
  end
  %
  disp('  test I/C/dV functions with vectorized input');
  if(n.ifc_simu && ~n.is_vsrc)
    n.I(rand(3,10));
    n.C(rand(3,10));
    n.dV(rand(3,10));
  end
  %
  disp('  test v for voltage source'); 
  if(v.is_vsrc), v.V(10); end

  disp('-----------------------------------------------------------------------');
  disp('Test verification interface');
  if(iro.ifc_verify)
    s = coho_shape(repmat([0,0.1],N,1),[]); 
    [c,err]=iro.I_ldi(s); 
    [c,err]=iro.dV_ldi(s);
    [c,err]=iro.dV_ldi(s,3);
  end

  disp('-----------------------------------------------------------------------');
  disp('Test small-signal-analysis interface');
  if(latch.ifc_ssan)
    j = latch.Jac(rand(latch.nodeNum,1));
    didv = latch.dIdV(rand(latch.nodeNum,1));
    % doesn't support array inputs
    %j = latch.Jac(rand(latch.nodeNum,10));
    %didv = latch.dIdV(rand(latch.nodeNum,10));
  end

  disp('-----------------------------------------------------------------------');
  disp('Test utilities');
  disp('  test is_leaf');
  assert(n.is_leaf); 
  assert(~latch.is_leaf);
 
  disp('  test nodeNum & elemNum');
  assert(latch.nodeNum==11&&latch.elemNum==5);

  disp('  test find_node_name');
  names = latch.find_node_name;
  for i=1:length(names)
    disp(names{i});
  end

  disp('  test find_port_index');
  i = latch.find_port_index(latch.d); 
  disp(latch.find_node_name(i));

  disp('    ');
  disp('  test print_circuit_tree'); 
  latch.print_circuit_tree;

  disp('    ');
  disp('  test print_nodes'); 
  latch.print_nodes; 

  disp('    ');
  disp('  test print_status'); 
  latch.print_status; 
end

function test_flat_vec
  disp('-----------------------------------------------------------------------');
  disp('Check the change of internal structures'); 
  clear classes
  addpath('./rro_v2'); 
  test_flat_vec_struct; 
  rmpath('./rro_v2');

  disp('-----------------------------------------------------------------------');
  disp('Check the results are similar'); 
  clear classes
  addpath('./rro_v2'); 
  test_flat_vec_result; 
  rmpath('./rro_v2');

  disp('-----------------------------------------------------------------------');
  disp('Compare the performance'); 
  disp('  Usually, the original circuit is the slowest, flattened is faster, flattened and vectorized is much faster');
  clear classes
  addpath('./rro_v1'); 
  test_flat_vec_perf; 
  rmpath('./rro_v1');

  clear classes
  disp(' However, vectorized but not flattened circuit could even be a little slower than the original circuit due to extra function calles');
  addpath('./rro_v2'); 
  test_flat_vec_perf; 
  rmpath('./rro_v2');
end

function test_flat_vec_struct
  rro = rambusOsc('rro',4,[1,1]*1e-5); 
  disp('The original circuit tree is');
  rro.print_circuit_tree('oo  '); 
  disp('The circuit tree after flatten is');
  rro.flatten.print_circuit_tree('ff  ');
  disp('The circuit tree after flatten and vectorize is');
  rro.vectorize.print_circuit_tree('vv  ');
  disp('The status after flatten and vectorize is');
  rro.vectorize.print_status;

  disp('Flatten and vectorize will change the circuit internal structure.'); 
  disp('Matlab copy is shallow copy. So to keep a copy, you may have to re-create the same circuit');
  % NOTE: Matlab's copy is shallow copy. So rro and rro1 shares the same rro_stage.  
  % rro1.vectorize will also change rro's elements' status.
  % So to keep a copy, just re-create the same circuit again
  rro = rambusOsc('rro1',2,[1,1]*1e-5);
  rro1 = copy(rro); rro1.vectorize; 
  rro.print_status;
end

function test_flat_vec_result
  % result should be similar, may not be the same because of different math op order.
  rro = rambusOsc('rro',4,[1,1]*1e-5);
  rro2 = copy(rro); rro2.flatten.flatten.vectorize;
  V = rand(4*2,10000);
  i1 = rro.I(V); i2 = rro2.I(V);
  c1 = rro.C(V); c2 = rro2.C(V);
  assert(max(abs(i1(:)-i2(:)))<10*eps);
  assert(max(abs(c1(:)-c2(:)))<10*eps);
end

% NOTE: waveform maybe different because of the difference in I() and other integration error
function test_flat_vec_perf
  do_opt = false; 
  rro1 = rambusOsc('rro1',4,[1,1]*1e-5);
  rro2 = rambusOsc('rro1',4,[1,1]*1e-5);
  rro3 = rambusOsc('rro1',4,[1,1]*1e-5);
  rro4 = rambusOsc('rro1',4,[1,1]*1e-5);
  rro2.flatten;
  rro3.vectorize;
  rro4.flatten.vectorize;
  s1 = testbench(rro1,{},{},do_opt);
  s2 = testbench(rro2,{},{},do_opt);
  s3 = testbench(rro3,{},{},do_opt);
  s4 = testbench(rro4,{},{},do_opt);

  v0 = rand(8,1)*1e-3;
  T = (0:0.001:1)*1e-9; 

  disp('Original circuit');
  tic
  [t1,v1,dv1]=s1.simulate(T,v0);
  toc
  figure; hold on; plot(t1,v1(:,1),'r');
  
  disp('Flattened circuit');
  tic
  [t2,v2,dv2]=s2.simulate(T,v0);
  toc
  plot(t2,v2(:,1),'g');
  

  disp('Vectorized circuit');
  tic
  [t3,v3,dv3]=s3.simulate(T,v0);
  toc
  plot(t3,v3(:,1),'b');

  disp('Flattened and vectorized circuit');
  tic
  [t4,v4,dv4]=s4.simulate(T,v0);
  toc
  plot(t4,v4(:,1),'k');
end

function test_src
  % isrc
  disp('Circuit class does support current source');
  isrc = test_isrc('test_isrc',1e-3,1e-9); 
  isrct = testbench(isrc,{},{});
  [t,v] = isrct.simulate([0,1e-5]);
  figure; plot(t,v);

  % vsrc
  disp('Circuit class does not support voltage source');
  %rc = test_vrc('rc',3,1e3,1e-9);
  %rct = testbench(rc,{},{}); 
  %[t,v] = rct.simulate([0,1e-5]); 
  %plot(t,v);

  disp('We have to make the node as input of circuits and apply the voltage current in testbench'); 
  rc = test_rc('rc',1e3,1e-9);
  rct = testbench(rc,{},{});
  v0 = 1; 
  [t,v] = rct.simulate([0,1e-5],v0); 
  figure; plot(t,v);

  % NOTE: why it's slow? 
  % RRO: 1000 points, 1 leaf, 10 seconds
  % NAND: 18K points, 3 leaves, 340 seconds
  disp('Another example: NAND with vdd/gnd specified in testbench. This may take more than 5 mins...');
  vnn = 1.8;
  vdd = vsrc('vdd',vnn); gnd = vsrc('gnd',0); 
  i1 = vsin('i1',1e9,0,1,vnn/2);
  i2 = vpulse('i2',[0,vnn],[1,4,1,4]*1e-9);
  nd = test_nand('nand',1e-5); 
  ndt = testbench(nd,{nd.vdd,nd.gnd,nd.i1,nd.i2},{vdd,gnd,i1,i2});
  [t,v] = ndt.simulate([0,2e-8]);
  figure; plot(t,v(:,3:end));
 
  % NOTE: this is slower because of less vectorization.  (nmos+snmos+spmos, vs. nmos+pmos)
  disp('Standard NAND. This may take more than 5 mins');
  nd2 = NAND('nand2',1e-5); 
  ndt2 = testbench(nd2,{nd2.i1,nd2.i2},{i1,i2}); 
  [t2,v2] = ndt2.simulate([0,2e-8]);
  figure; plot(t2,v2); 
  assert(max(max(abs(v2-v(:,(3:end)))))==0);
end
