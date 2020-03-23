function plotCap(data)

Vd = data.Vd;
Vg = data.Vg;
Vs = data.Vs;
Vb = data.Vb;

Cdb = data.Cdb;
Csb = data.Csb;
Cgb = data.Cgb;

figure
p1 = subplot(3,1,1);
hold on
p2 = subplot(3,1,2);
hold on
p3 = subplot(3,1,3);
hold on

for i = 1:11
    Vdi = Vd(1+(i-1)*101:i*101);
    Vsi = Vs(1+(i-1)*101:i*101);
    Vgi = Vg(1+(i-1)*101:i*101);
    Vbi = Vb(1+(i-1)*101:i*101);
    Cdbi = Cdb(1+(i-1)*101:i*101);
    Csbi = Csb(1+(i-1)*101:i*101);
    Cgbi = Cgb(1+(i-1)*101:i*101);
    
    plot(p1,Vdi-Vbi,Cdbi)
    plot(p2,Vsi-Vbi,Csbi)
    plot(p3,Vgi-Vbi,Cgbi);

end

