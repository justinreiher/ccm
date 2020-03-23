Simulation PTM Models


.include 45nm_HP.pm

Vd vd gnd DC 0
Vg vg gnd DC 0
Vs vs gnd DC 1.0
Vb vb gnd DC 1.0

*Mname NodeDrain NodeGate NodeSource NodeBody
M1 vd vg vs vb pmos l = 45n w = 450n

* Generates data for Id-Vd curves at different Vgs
.dc vd 0 1.0 0.01 vg 0 1.0 0.1

.control
save @m1[cgs] @m1[csg] @m1[cgd] @m1[cdg] @m1[cgb] @m1[cbg] @m1[csb] @m1[cbs] @m1[cdb] @m1[cbd] @m1[id] all
run
wrdata cpmosVdsSweep.txt (abs(@m1[cgs])+abs(@m1[csg])) (abs(@m1[cgd])+abs(@m1[cdg])) (abs(@m1[cgb])+abs(@m1[cbg])) (abs(@m1[csb])+abs(@m1[cbs])) (abs(@m1[cdb])+abs(@m1[cbd])) @m1[id]
wrdata pmos45nmHP_VdsSweep.txt -i(vd) v(vd) v(vg) v(vs) v(vb)
plot -i(vd)
*show m1

.endc


.end
