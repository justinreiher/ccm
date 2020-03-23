function params = fit_ekv(startParams,sweepY,sweepG,sweepX,sweeps)
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

% remove anywhere where Vx-Vb == 0, derivative is not defined at that point
ind = (VxS-VbS) == 0;
VyS = VyS(~ind);
VgS = VgS(~ind);
VxS = VxS(~ind);
VbS = VbS(~ind);
IdS = IdS(~ind);


biasVoltage = [VyD,VyG;     %Drain  (Vy)
                VgD,VgG;    %Gate   (Vg)
                VxD,VxG;    %Source (Vx)
                VbD,VbG];   %Body   (Vb)
            
biasVoltageBodyFit = [VyS;VgS;VxS;VbS];

Id_ref = [IdD,IdG];

fitParamsPrev = zeros(size(fitParams));
weight = computeWeights(biasVoltage,lengthRuns);

paramPrev = populateParamStruct(fitParamsPrev,startParams);
param = populateParamStruct(fitParams,startParams);

Id = ekv_Ids(param,biasVoltage);
gradParams_Ids = gradParams_ekvIds(paramPrev,biasVoltage);
gradErrPrev = sum(gradParams_Ids*diag(2*weight.*(Id-Id_ref)),2);

tol = 1e-12;


params = populateParamStruct(fitParams,startParams);
err = (norm(Id-Id_ref)/norm(Id_ref));
errPrev = 0;
while( abs(err-errPrev) > tol)
    
    fprintf('Relative Error: %d\n',err);

    gradParams_Ids = gradParams_ekvIds(params,biasVoltage);
    gradErr = sum(gradParams_Ids*diag(2*weight.*(Id-Id_ref)),2);
    
    %wolfe condition
    gamma = abs((fitParams - fitParamsPrev)'*(gradErr - gradErrPrev))/(norm(gradErr-gradErrPrev))^2;
    
    fitParamsPrev = fitParams;
    gradErrPrev = gradErr;
    fitParams = fitParams - gamma*gradErr;
    
    params = populateParamStruct(fitParams,startParams);
    Id = ekv_Ids(params,biasVoltage);
    errPrev = err;
    err = norm(Id-Id_ref)/(norm(Id_ref));

end

fitParams(5) = 0.2;
fitParams(6) = 0.5;
weight = eye(length(IdS));
Id = ekv_Ids(param,biasVoltageBodyFit);
gradParams_Ids = gradParams_ekvIdsBody(paramPrev,biasVoltageBodyFit);
gradErrPrev = sum(gradParams_Ids*diag(2*weight.*(Id-IdS)),2);
params = populateParamStruct(fitParams,startParams);
errPrev = 0;
while( abs(err-errPrev) > tol)
    
    fprintf('Relative Error: %d\n',err);

    gradParams_Ids = gradParams_ekvIdsBody(params,biasVoltageBodyFit);
    gradErr = sum(gradParams_Ids*diag(2*weight.*(Id-IdS)),2);
    
    %wolfe condition
    gamma = abs((fitParams - fitParamsPrev)'*(gradErr - gradErrPrev))/(norm(gradErr-gradErrPrev))^2;
    
    fitParamsPrev = fitParams;
    gradErrPrev = gradErr;
    fitParams = fitParams - gamma*gradErr;
    
    params = populateParamStruct(fitParams,startParams);
    Id = ekv_Ids(params,biasVoltageBodyFit);
    errPrev = err;
    err = norm(Id-IdS)/(norm(IdS));

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

% if type == 'p'
%     addedWeight = find(V(3,:) > 0.5);
%     weights(addedWeight) = 1000;
% end
% 
% if type == 'n'
%     addedWeight = find(V(1,:) > 0.5);
%     weights(addedWeight) = 1000;
% end

end
