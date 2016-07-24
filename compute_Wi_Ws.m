function [Wi, Ws] = compute_Wi_Ws(I0, VX, VY, h2d, TP, params)
hsz = TP.HSZ;
S   = TP.S;
Wi  = cell(2*params.temporalR+1,1);
i   = -params.nfore:params.nback;

for j = 1:length(i)
    vx = VX(:, j);
    vy = VY(:, j);
    J  = TP.NSY(:, j);

    Wi{j} = approxL1inv(I0, J, vx, vy, S, h2d, params.epsilon, hsz);
end

Ws = approxL1invWS(I0, params.epsilon, hsz);
end