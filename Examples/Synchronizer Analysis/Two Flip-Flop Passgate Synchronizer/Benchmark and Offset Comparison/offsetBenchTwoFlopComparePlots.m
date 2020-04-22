%%% Strings
% this script intializes the functions to plot and analyze the differences
% between a 2 flip-flop passgate design with all transistors of the same
% width and one with an offset in the coupling inverter between stages 2
% and 3 with a 2:1 ratio. This line may be commented out after the first
% run through
componentGainAnalysis2Flop

syncRef = TWO_FLOP_PG_SYNC('ref',1,1);

inv3str = {'$q_0$','$\overline{q_0}$','inv3'}; %coupling inverter

pg2str  = {'$\overline{q_0}$','$x_1$','$pg_2$'}; %coupling pg
pg2Nstr = {'$\overline{q_0}$','$x_1$','$pg_2$ - NMOS'};
pg2Pstr = {'$\overline{q_0}$','$x_1$','$pg_2$ - PMOS'};

inv4str = {'$x_1$','$q_1$','inv4'};            %fwd
inv5str = {'$q_1$','$y_1$','inv5'};            %bwk

pg3str  = {'$y_1$','$x_1$','$pg_3$'};            %feedback pg
pg3Nstr = {'$y_1$','$x_1$','$pg_3$ - NMOS'};
pg3Pstr = {'$y_1$','$x_1$','$pg_3$ - PMOS'};


startTime = 6e-10;
endTime   = 1e-9;

%%
% coupling Inverter 3
inv3In = syncRef.cctPath.PGFF_0.PGL_0.q;
inv3Out = syncRef.cctPath.PGFF_0.PGL_0.qbar;

[inv3O,t_inv3O] = invGOffset(inv3In,inv3Out,[startTime,endTime]);
[inv3B,tinv3B] = invGBench(inv3In,inv3Out,[startTime,endTime]);
componentCompare(inv3In,inv3Out,tinv3B,inv3B,t_inv3O,inv3O,inv3str);

%%
% passgate 2 - coupling gate
pg2In  = syncRef.cctPath.PGFF_0.PGL_0.qbar;
pg2Out = syncRef.cctPath.PGFF_0.PGL_1.x0;

[pg2O,tpg2O] = invGOffset(pg2In,pg2Out,[startTime,endTime]);
[pg2B,tpg2B] = invGBench(pg2In,pg2Out,[startTime,endTime]);

[pg2ON,tpg2ON] = txGOffset(pg2In,pg2Out,[startTime,endTime],'nmos');
[pg2BN,tpg2BN] = txGBench(pg2In,pg2Out,[startTime,endTime],'nmos');

[pg2OP,tpg2OP] = txGOffset(pg2In,pg2Out,[startTime,endTime],'pmos');
[pg2BP,tpg2BP] = txGBench(pg2In,pg2Out,[startTime,endTime],'pmos');


componentCompare(pg2In,pg2Out,tpg2B,pg2B,tpg2O,pg2O,pg2str);
componentCompare(pg2In,pg2Out,tpg2BN,pg2BN,tpg2ON,pg2ON,pg2Nstr);
componentCompare(pg2In,pg2Out,tpg2BP,pg2BP,tpg2OP,pg2OP,pg2Pstr);

%%
% inverter 4 - Forward CC
inv4In  = syncRef.cctPath.PGFF_0.PGL_1.x0;
inv4Out = syncRef.cctPath.PGFF_0.q;

[inv4O,tinv4O] = invGOffset(inv4In,inv4Out,[startTime,endTime]);
[inv4B,tinv4B] = invGBench(inv4In,inv4Out,[startTime,endTime]);

componentCompare(inv4In,inv4Out,tinv4B,inv4B,tinv4O,inv4O,inv4str);

%%
% inverter 5 - Backward CC
inv5In  = syncRef.cctPath.PGFF_0.q;
inv5Out = syncRef.cctPath.PGFF_0.PGL_1.y0;

[inv5O,tinv5O] = invGOffset(inv5In,inv5Out,[startTime,endTime]);
[inv5B,tinv5B] = invGBench(inv5In,inv5Out,[startTime,endTime]);

componentCompare(19,22,tinv5B,inv5B,tinv5O,inv5O,inv5str);

%%
% passgate 3 - feedback PG
pg3In  = syncRef.cctPath.PGFF_0.PGL_1.y0;
pg3Out = syncRef.cctPath.PGFF_0.PGL_1.x0;

[pg3O,tpg3O] = invGOffset(pg3In,pg3Out,[startTime,endTime]);
[pg3B,tpg3B] = invGBench(pg3In,pg3Out,[startTime,endTime]);

[pg3ON,tpg3ON] = txGOffset(pg3In,pg3Out,[startTime,endTime],'nmos');
[pg3BN,tpg3BN] = txGBench(pg3In,pg3Out,[startTime,endTime],'nmos');

[pg3OP,tpg3OP] = txGOffset(pg3In,pg3Out,[startTime,endTime],'pmos');
[pg3BP,tpg3BP] = txGBench(pg3In,pg3Out,[startTime,endTime],'pmos');

componentCompare(pg3In,pg3Out,tpg3B,pg3B,tpg3O,pg3O,pg3str);
componentCompare(pg3In,pg3Out,tpg3BN,pg3BN,tpg3ON,pg3ON,pg3Nstr);
componentCompare(pg3In,pg3Out,tpg3BP,pg3BP,tpg3OP,pg3OP,pg3Pstr);


