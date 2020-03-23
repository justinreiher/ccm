function Cdevs = ekvModelCap(V,capParams)
%ekvModelCap computes the capacitance matrix for each device based on the
%voltage on the terminals of the device

i_f = capParams.i_f;
i_r = capParams.i_r;
n   = capParams.n;

q_f = (sqrt(1 + 4*i_f) - 1)/2;
q_r = (sqrt(1 + 4*i_r) - 1)/2;

junctionCap = capParams.junctionCap;
cox = capParams.cap;

[terminalNum,numDevices] = size(V);

Vdb = V(1,:) - V(4,:);
Vsb = V(3,:) - V(4,:);

c_s = (1/3)*(q_f.*(2*q_f+4*q_r + 3))./(q_f + q_r +1).^2;
c_d = (1/3)*(q_r.*(2*q_r+4*q_f + 3))./(q_r + q_f +1).^2;
Cgs = cox.*c_s;
Cgd = cox.*c_d;
Cgb = ((n-1)./n).*(cox -Cgs-Cgd);



Cbd = junctionDrainCap(Vdb,junctionCap); %+ (n-1).*Cgs;
Cbs = junctionDrainCap(Vsb,junctionCap); %+ (n-1).*Cgd;


Cdevs = repmat(Vdb*0,16,1);

%drain caps
Cdevs(5,:)  = Cgd;
Cdevs(13,:) = Cbd;
%gate caps
Cdevs(2,:)  = Cgd;
Cdevs(10,:) = Cgs;
Cdevs(14,:) = Cgb;
%source caps
Cdevs(7,:)  = Cgs;
Cdevs(15,:) = Cbs;
%body caps
Cdevs(4,:)  = Cbd;
Cdevs(8,:)  = Cgb;
Cdevs(12,:) = Cbs;

%Cdev = [0, Cgd,0.125,0.125;
%       0.25,0,0.25,0.25;
%    0.125,0.25,0,0.125;
%    0.125,0.25,0.125,0];

%Cdevs = repmat(Cdev,1,numDevices);
Cdevs = reshape(Cdevs,terminalNum,terminalNum,[]);

end

