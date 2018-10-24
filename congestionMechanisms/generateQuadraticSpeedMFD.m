function congestion = generateQuadraticSpeedMFD(C,L,vf)
% This function creates a stucture with various features of a quadratic
% MFD.
% With a quadratic speed MFD, C=ncr*vcr/L=ncr vf (1-ncr/(3ncr))^2 /L
% =ncr vf (2/3)^2 /L => ncr= 9/4 C L /vf

% Last modified by Raphael Lamotte, on October 24, 2018.
congestion.mechanism='MFD';
congestion.C=C; % the maximum trip completion rate possible in steady state.
congestion.vf=vf; % Free-flow speed
congestion.ncr=9/4*C*mean(L)/vf; % accumulation corresponding to C in steady state.
congestion.speed_fct=@(x)vf*max(sign(1-x/(3*congestion.ncr)).*(1-x/(3*congestion.ncr)).^2,0.05); % speed-MFD
end

