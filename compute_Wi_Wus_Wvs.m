function [Wi, Wus, Wvs] = compute_Wi_Wus_Wvs(I0, J, vx, vy, h2d, TP, params)
    hsz  = TP.HSZ;
    S    = TP.S;

    %% Wi
    Wi   = approxL1inv(I0, J, vx, vy, S, h2d, params.epsilon, hsz);
    %% Wus
    Wus  = approxL1invWS(vx, params.epsilon, hsz);
    %% Wvs
    Wvs  = approxL1invWS(vy, params.epsilon, hsz);
end