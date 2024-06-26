function state = dynamics(~, S, Xddotf, Yddotf, Zddotf, Af, const)

% Define variables to be used in computation
thetadot = const.thetadot;
Re = const.Re;
Ar = const.Ar;
m = const.m;
rho0 = const.rho0;
H = const.H;
r0 = const.r0;


% Evaluate Acceleration terms
aX = Xddotf(Ar, S(9), H, S(8), Re, S(1), S(4), S(2), S(5), S(3), S(6),...
    m, S(7), r0, rho0, thetadot);
aY = Yddotf(Ar, S(9), H, S(8), Re, S(1), S(4), S(2), S(5), S(3), S(6),...
    m, S(7), r0, rho0, thetadot);
aZ = Zddotf(Ar, S(9), H, S(8), Re, S(1), S(4), S(2), S(5), S(3), S(6),...
    m, S(7), r0, rho0, thetadot);

% Evaluate the A matrix
Ac = Af(Ar, S(9), H, S(8), Re, S(1), S(4), S(2), S(5), S(3), S(6),...
    m, S(7), r0, rho0, thetadot);

% Reshape Phi for easier use
phi = reshape(S(19:end), 18, 18);

% define final state 
state = [S(4); S(5); S(6); aX; aY; aZ; 0; 0; 0; 0; 0;...
         0; 0; 0; 0; 0; 0; 0; reshape(Ac*phi, 324, 1)];

end