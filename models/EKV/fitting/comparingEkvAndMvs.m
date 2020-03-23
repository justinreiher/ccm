Vd = 0:0.01:1;
Vs = zeros(size(Vd));
VbN = zeros(size(Vd));
VbP = ones(size(Vd));
Vg = 0:0.1:1;
VdP = ones(size(Vd));
VgSweep = 0:0.01:1;


pmosMVS = mvs_model_params_45nmHP('p');
nmosMVS = mvs_model_params_45nmHP('n');

pmosEKV = ekv_model_params_45nmHP('p');
nmosEKV = ekv_model_params_45nmHP('n');


figure
pmos = subplot(1,1,1);
title('pmos V_{ds}')
hold on

figure
nmos = subplot(1,1,1);
title('nmos V_{ds}')
hold on

figure
nmosVG = subplot(1,1,1);
title('nmos V_{gs}')
hold on

figure
pmosVG = subplot(1,1,1);
title('pmos V_{gs}')
hold on

for i = 1:length(Vg)
    plot(pmos,Vd-Vs,mvs_Id(pmosMVS,[Vd;Vg(i)*ones(size(Vd));Vs;VbP]),'r')
    plot(pmos,Vd-Vs,ekv_Ids(pmosEKV,[Vd;Vg(i)*ones(size(Vd));Vs;VbP]))
    plot(pmos,Vs-Vd,mvs_Id(pmosMVS,[Vs;Vg(i)*ones(size(Vd));Vd;VbP]),'r')
    plot(pmos,Vs-Vd,ekv_Ids(pmosEKV,[Vs;Vg(i)*ones(size(Vd));Vd;VbP]))
    

    
    plot(nmos,Vd-Vs,mvs_Id(nmosMVS,[Vd;Vg(i)*ones(size(Vd));Vs;VbN]),'r')
    plot(nmos,Vd-Vs,ekv_Ids(nmosEKV,[Vd;Vg(i)*ones(size(Vd));Vs;VbN]))
    plot(nmos,Vs-Vd,mvs_Id(nmosMVS,[Vs;Vg(i)*ones(size(Vd));Vd;VbN]),'r')
    plot(nmos,Vs-Vd,ekv_Ids(nmosEKV,[Vs;Vg(i)*ones(size(Vd));Vd;VbN]))
    
    
end

    plot(pmosVG,VgSweep-Vs,mvs_Id(pmosMVS,[VdP;VgSweep;Vs;VbP]),'r')
    plot(pmosVG,VgSweep-Vs,ekv_Ids(pmosEKV,[VdP;VgSweep;Vs;VbP]))
    plot(pmosVG,VgSweep-VdP,mvs_Id(pmosMVS,[Vs;VgSweep;VdP;VbP]),'r')
    plot(pmosVG,VgSweep-VdP,ekv_Ids(pmosEKV,[Vs;VgSweep;VdP;VbP]))

    plot(nmosVG,VgSweep-Vs,mvs_Id(nmosMVS,[VdP;VgSweep;Vs;VbN]),'r')
    plot(nmosVG,VgSweep-Vs,ekv_Ids(nmosEKV,[VdP;VgSweep;Vs;VbN]))
    plot(nmosVG,VgSweep-VdP,mvs_Id(nmosMVS,[Vs;VgSweep;VdP;VbN]),'r')
    plot(nmosVG,VgSweep-VdP,ekv_Ids(nmosEKV,[Vs;VgSweep;VdP;VbN]))