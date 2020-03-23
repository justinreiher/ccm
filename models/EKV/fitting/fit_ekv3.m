function params = fit_ekv(startParams,sweepY,sweepG,sweepX,sweeps,type)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here



%%%%%%%%%%%%%%%%%%%
%          Vy (Drain)
%         _|
%      | |
% Vg --| |-| Vb (not considered in this model)
%      | |_ 
%          |
%          Vx (Source)
%%%%%%%%%%%%%%%%%%%

%fine grain sweeping
fineStart = sweeps(2);
fineStep = sweeps(3);
fineEnd = sweeps(4);

%coarse step sweeping
coarseStart = sweeps(5);
coarseStep  = sweeps(6);
coarseEnd   = sweeps(7);


lengthRuns = length(fineStart:fineStep:fineEnd);
numSteps   = length(coarseStart:coarseStep:coarseEnd);

fitParams = [startParams.Io;
             startParams.alpha;
             startParams.beta;
             startParams.Vth;
             startParams.gamma;
             startParams.phi];

if type == 'p'
    fitParams = -fitParams;
    fitParams(5) = -fitParams(5);
    fitParams(6) = -fitParams(6);
end

VyD = sweepY(:,1)';
VgD = sweepY(:,2)';
VxD = sweepY(:,3)';
VbD = sweepY(:,4)';
IdD = sweepY(:,5)';

VyG = sweepG(:,1)';
VgG = sweepG(:,2)';
VxG = sweepG(:,3)';
VbG = sweepG(:,4)';
IdG = sweepG(:,5)';

VyS = sweepX(:,1)';
VgS = sweepX(:,2)';
VxS = sweepX(:,3)';
VbS = sweepX(:,4)';
IdS = sweepX(:,5)';


biasVoltage = [VyD,VyG,VyS;     %Drain  (Vy)
                VgD,VgG,VgS;    %Gate   (Vg)
                VxD,VxG,VxS;    %Source (Vx)
                VbD,VbG,VbS];   %Body   (Vb)
            

Id_ref = [IdD,IdG,IdS];

fitParamsPrev = zeros(size(fitParams));
% phi cannot start at zero
fitParamsPrev(6) = 0.7;

tol = 1e-12;

paramPrev = populateParamStruct(fitParamsPrev,startParams);
param = populateParamStruct(fitParams,startParams);

weight = computeWeights(biasVoltage,lengthRuns,type);
Id = ekv_Ids(param,biasVoltage);
gradParams_Ids = gradParams_ekvIdsBody(paramPrev,biasVoltage);
gradErrPrev = sum(gradParams_Ids*diag(2*weight.*(Id-Id_ref)),2);
params = populateParamStruct(fitParams,startParams);
err = (norm(Id-Id_ref)/norm(Id_ref));
errPrev = 0;
while( abs(err-errPrev) > tol)
    
    fprintf('Relative Error: %d\n',err);

    gradParams_Ids = gradParams_ekvIdsBody(params,biasVoltage);
    gradErr = sum(gradParams_Ids*diag(2*weight.*(Id-Id_ref)),2);
    
    %wolfe condition
    gamma = abs((fitParams - fitParamsPrev)'*(gradErr - gradErrPrev))/(norm(gradErr-gradErrPrev))^2;
    
    fitParamsPrev = fitParams;
    gradErrPrev = gradErr;
    fitParams = fitParams - gamma*gradErr;
    
    if fitParams(6) < 0
        fitParams(6) = fitParamsPrev(6)/2;
    end
    
    params = populateParamStruct(fitParams,startParams);
    Id = ekv_Ids(params,biasVoltage);
    errPrev = err;
    err = norm(Id-Id_ref)/(norm(Id_ref));

end


end


function params = populateParamStruct(input,ref)

    params = ref;
    params.Io = input(1);
    params.alpha = input(2);
    params.beta = input(3);
    params.Vth = input(4);
    params.gamma = input(5);
    params.phi = input(6);
end

function weights = computeWeights(V,lengthRuns,type)

%%%%%%%%%%
% V(1,:) = all of Vd
% V(2,:) = all of Vg
% V(3,:) = all of Vs
% V(4,:) = all of Vb
%%%%%%%%

weights = ones(1,length(V));

 if type == 'p'
     addedWeight = find(V(3,:) > 0.5);
     weights(addedWeight) = 1000;
 end
 
 if type == 'n'
     addedWeight = find(V(1,:) > 0.5);
     weights(addedWeight) = 1000;
 end

end
