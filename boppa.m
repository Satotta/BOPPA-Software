function [RMS, x0_hat] = boppa(rv0, site, tspan, x0_bar, P0_bar, R, obs_data, n, const)
% The following user defined function uses the Batch Processor algorithm in
% the orbit-determination procedure. Orbit dynamics are based on the
% two-body equation of motion with included effects from Earth's oblateness
% (J2) and atmospheric drag. The function will run the algorithm n times,
% after which convergence is expected. The function can handle observations
% from up to three different sites whose coordinates will be specified by
% the user. The inputs and outputs are described below.


%-------------------------------------------------------------------------%
%                              Inputs Guide                               %
%-------------------------------------------------------------------------%

% --> rv0: a 6x1 column array containing the spacecraft initial position [m]
% (rows 1:3) and velocity [m/s] (rows 4:6) in the ECI frame

% --> site: a 3x3 matrix containing the position [m] of up to three sites;
% each column corresponds to a particular site with the X, Y, and Z
% positions within the ECEF frame being specified by rows 1 through 3
% respectively

% --> tspan: a time array [s] of variable length specified by the user in
% the form t0:step:tf

% --> x0_bar: an 18x1 column vector holding the apriori state fo the
% spacecraft

% --> P0_bar: an 18x18 covariance matrix

% --> R: an 2x2 matrix that reflects the noise on the data

% --> obs_data: a matrix of variable length containing the time stamp for
% each observation in column 1, the site name in column 2, and the range
% and range rate observation in columns 3 and 4 respectively

% --> n: number of batch iterations (3 recommended)

% --> const: a structure array containing the dynamic and reference
% constants to be used in calculation [units vary]


%-------------------------------------------------------------------------%
%                            Outputs Guide                                %
%-------------------------------------------------------------------------%

% --> RMS: a 2xn matrix holding the range and range rate residual RMS value
% for each batch iteration. The first row corresponds to the range residual
% RMS values and the second row corresponds to the range rate residual RMS
% values

% --> x0_hat: an 18xn+1 matrix holding the state updates mapped back to the
% initial epoch. The nth column corresponds to the nth batch iteration and 
% the last column is the total state change from the initial state.

% --> Residual Plots: a range and range rate residual plot for each
% iteration will be outputted

% --> Position Deviation Error Ellipsoid: a 3D position deivation ellipsoid
% will be outputted

% --> Velocity Deviation Error Ellipsoid: a 3D velocity deivation ellipsoid
% will be outputted




%-------------------------------------------------------------------------%
%     Symbolic Determination of Important Batch Processor Parameters      %
%-------------------------------------------------------------------------%


% Define the symbolic variables to be used in subsequent parameter
% determination
syms t X Y Z Xdot Ydot Zdot mu J2 Cd Xs Ys Zs Xs1 Ys1 Zs1 ...
    Xs2 Ys2 Zs2 Xs3 Ys3 Zs3 thetadot Re rho0 H r0 Ar m


% Define the equations for range and range rate for an arbitrary site
theta = thetadot*t;
c = cos(theta);
s = sin(theta);

rho = simplify(sqrt((X - (Xs*c - Ys*s))^2 + ...
                    (Y - (Xs*s + Ys*c))^2 + ...
                    (Z - Zs)^2));

rhodot = simplify(1/rho*(Xdot*(X - Xs*c + Ys*s) + ...
                               Ydot*(Y - Ys*c - Xs*s) + ...
                               Zdot*(Z - Zs) + thetadot*...
                               ((X*Xs + Y*Ys)*s + ...
                               (X*Ys - Y*Xs)*c)));


% Develop the symbolic range and range rate matrix for an arbitrary site
O = [rho; rhodot];


% Determine the symbolic state for each site
S_101 = [X, Y, Z, Xdot, Ydot, Zdot, mu, J2, Cd, ...
    Xs, Ys, Zs, Xs2, Ys2, Zs2, Xs3, Ys3, Zs3];

S_337 = [X, Y, Z, Xdot, Ydot, Zdot, mu, J2, Cd, ...
    Xs1, Ys1, Zs1, Xs, Ys, Zs, Xs3, Ys3, Zs3];

S_394 = [X, Y, Z, Xdot, Ydot, Zdot, mu, J2, Cd, ...
    Xs1, Ys1, Zs1, Xs2, Ys2, Zs2, Xs, Ys, Zs];


% Determine the Ht matrix for each site
Ht_101 = simplify(jacobian(O, S_101));
Ht_337 = simplify(jacobian(O, S_337));
Ht_394 = simplify(jacobian(O, S_394));


% Define symbolic equations to represent the two-body equation of motion
% and J2 dynamics
r = sqrt(X^2 + Y^2 + Z^2);
Uprime = mu/r*(1 - J2*(Re/r)^2*(3/2*(Z/r)^2 - 1/2));
aJ2_X = diff(Uprime, X);
aJ2_Y = diff(Uprime, Y);
aJ2_Z = diff(Uprime, Z);


% Define a symbolic equation for the drag dynamics
vA_v = [Xdot + thetadot*Y; Ydot - thetadot*X; Zdot];
vA = norm(vA_v);
rho_D = rho0*exp(-(r - r0)/H);
aD = -1/2*Cd*Ar/m*rho_D*vA.*vA_v;


% Use in conjunction to determine the total dynamical equations in
% component form
Xddot = simplify(aJ2_X + aD(1));
Yddot = simplify(aJ2_Y + aD(2));
Zddot = simplify(aJ2_Z + aD(3));


% Define the dynamics in state space form
F = [Xdot; Ydot; Zdot; Xddot; Yddot; Zddot; ...
    0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];


% Define the state vector
S_A = [X, Y, Z, Xdot, Ydot, Zdot, mu, J2, Cd, ...
    Xs1, Ys1, Zs1, Xs2, Ys2, Zs2, Xs3, Ys3, Zs3];


% Determine the A matrix via the jacobian
A = simplify(jacobian(F, S_A))



%-------------------------------------------------------------------------%
%                       Batch Processor Algorithm                         %
%-------------------------------------------------------------------------%


% Convert the symbolic expressions into matlab functions to be used in
% subsequent calculations
rhof = matlabFunction(rho);
rhodotf = matlabFunction(rhodot);
Xddotf = matlabFunction(Xddot);
Yddotf = matlabFunction(Yddot);
Zddotf = matlabFunction(Zddot);
Af = matlabFunction(A);
Ht_101f = matlabFunction(Ht_101);
Ht_337f = matlabFunction(Ht_337);
Ht_394f = matlabFunction(Ht_394);


% Form complete initial condition for use in ODE45
X0 = [rv0; const.mu; const.J2; const.Cd; reshape(site, 9, 1);...
      reshape(eye(18), 324, 1)];


% Define the initial normal and lambda matrices
L = P0_bar\eye(18);
N = L*x0_bar;
W = R\eye(2);


% Sort imported data for easier use
t_o = obs_data(:, 1);
site_num = obs_data(:, 2);
Y = [obs_data(:, 3)'; obs_data(:, 4)'];


% Accumulate current observation for each time
% Define variables to be used in loop
y = zeros(2, length(t_o), n);
rho_c = zeros(1, length(t_o));
rhodot_c = zeros(1, length(t_o));
x0_hat = zeros(18, n);
P0 = zeros(18, 18, n);
RMS = zeros(2, n);
count = 1;


% Define loop to iterate enough for convergence or to terminate after 2
% iterations
for i = 1:n
    
    % Call ODE45 to integrate the equations of motion and update the state
    % transition matrix
    options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    [~, Ss] = ode45(@(t, S) dynamics(t, S, Xddotf, Yddotf, Zddotf, Af, const), ...
        tspan, X0, options);
    
    % Reshape the output for ease of use
    S = reshape(Ss', 18, 19, length(tspan));
    % The output corresponds to n 18x19 matrices where each matrix
    % corresponds to a time step specified by tspan. The first column
    % corresponds to the spacecraft state and columns 2 through 19
    % correspond to the STM.

    % Pull out state and STM arrays/Matrices
    X = S(:, 1, :);
    Phi = S(:, 2:end, :);

    % Accumulate current observation for each time
    for j = 1:length(t_o)

        % Define the position of the current observation time within the
        % time span used for integration
        k = find(tspan == t_o(j));

        % Use the index k to define the position and velocity at the
        % current observation time
        Xi = X(1, 1, k);
        Yi = X(2, 1, k);
        Zi = X(3, 1, k);
        Xdoti = X(4, 1, k);
        Ydoti = X(5, 1, k);
        Zdoti = X(6, 1, k);
        
        % Define switch statement to determine corresponding site
        switch site_num(j)
    
            case 101
                % Define site coordinates and Htilde matrix
                Xsi = X(10, 1, k);
                Ysi = X(11, 1, k);
                Zsi = X(12, 1, k);
                Ht = Ht_101f(Xi, Xdoti, Xsi, Yi, Ydoti, Ysi,...
                    Zi, Zdoti, Zsi, t_o(j), const.thetadot);
    
            case 337
                % Define site coordinates and Htilde matrix
                Xsi = X(13, 1, k);
                Ysi = X(14, 1, k);
                Zsi = X(15, 1, k);
                Ht = Ht_337f(Xi, Xdoti, Xsi, Yi, Ydoti, Ysi,...
                    Zi, Zdoti, Zsi, t_o(j), const.thetadot);
    
            otherwise
                % Define site coordinates and Htilde matrix
                Xsi = X(16, 1, k);
                Ysi = X(17, 1, k);
                Zsi = X(18, 1, k);
                Ht = Ht_394f(Xi, Xdoti, Xsi, Yi, Ydoti, Ysi,...
                    Zi, Zdoti, Zsi, t_o(j), const.thetadot);
        end 
      
        % Define the H matrix via the Ht matrix and STM matrix for the
        % current observation time
        H = Ht*Phi(:, :, k);
        
        % Determine calculated range and range rate for the given site
        rho_c(j) = rhof(Xi, Xsi, Yi, Ysi, Zi, Zsi, t_o(j), const.thetadot);
        rhodot_c(j) = rhodotf(Xi, Xdoti, Xsi, Yi, Ydoti, Ysi,...
            Zi, Zdoti, Zsi, t_o(j), const.thetadot);

        % Calculate residuals
        y(:, j, i) = Y(:, j) - [rho_c(j); rhodot_c(j)];
        
        % Define the updated lambda vector
        L = L + H'*W*H;

        % Define the updated normal matrix
        N = N + H'*W*y(:, j, i);

    end

    % Find best estimate using batch
    x0_hat(:, i) = L\N;

    % Find new covariance matrix
    P0(:, :, i) = L\eye(18);

    % Update nominal state
    X0 = [X0(1:18) + x0_hat(:, i); X0(19:end)];

    % Update apriori
    x0_bar = x0_bar - x0_hat(:, i);

    % Reset lambda
    L = P0_bar\eye(18);

    % Determine new normal matrix based on reset lambda and new x0_bar
    N = L*x0_bar;
    
    % Computs the RMS of the residuals for the current iteration
    RMS(:, i) = [sqrt(sum(y(1, :, i).^2)/length(t_o)) ;
                 sqrt(sum(y(2, :, i).^2)/length(t_o))];



%-------------------------------------------------------------------------%
%                             Residuals Plot                              %
%-------------------------------------------------------------------------%

    % Plot the residuals
    resFig = figure(1);
    subplot(n, 2, count)
    plot(y(1, :, i), 'LineWidth', 1.5), grid, xlabel('Observation Number'), ...
        ylabel('Range Residual (m)'), title(['Range Res. Pass ', num2str(i)])
    subplot(n, 2, count + 1)
    plot(y(2, :, i), 'LineWidth', 1.5), grid, xlabel('Observation Number'), ...
        ylabel('Range Rate Residual (m/s)'), ...
        title(['Range Rate Res. Pass ', num2str(i)])
    
    % Update plot count
    count = count + 2;
  
end


% Organize the best estimate of the initial state update after every pass
x0_hat = [x0_hat, sum(x0_hat, 2)];



%-------------------------------------------------------------------------%
%                    Position Error Ellipsoid Plot                        %
%-------------------------------------------------------------------------%

% Plot the 3D Position ellipsoid from the covariance corresponding to the last
% iteration

% Pull out corresponding batch covariance
Ppos = P0(1:3, 1:3, n);

% Determine eigen values and vectors
[vecsPos, valsPos] = eig(Ppos);

% Generate the general ellipsoid points for plotting
[ellipPos_x, ellipPos_y, ellipPos_z] = ellipsoid(0, 0, 0, 1, 1, 1);

% Scale the ellipsoid points w.r.t eigenvalues
ellipVel = [ellipPos_x(:), ellipPos_y(:), ellipPos_z(:)]*sqrt(valsPos)*vecsPos';
ellipPos_x(:) = ellipVel(:, 1);
ellipPos_y(:) = ellipVel(:, 2);
ellipPos_z(:) = ellipVel(:, 3);

% Plot the ellipsoid points
ellipPosFig = figure;
hold on
subplot(1, 2, 1)
surf(ellipPos_x, ellipPos_y, ellipPos_z), xlabel('X-Pos Deviation (m)'), ...
    ylabel('Y-Pos Deviation (m)'), zlabel('Z-Pos Deviation (m)'), ...
    title('Position Deviation Error Ellipsoid')



%-------------------------------------------------------------------------%
%                    Velocity Error Ellipsoid Plot                        %
%-------------------------------------------------------------------------%

% Plot the 3D Velocity ellipsoid from the covariance corresponding to the last
% iteration

% Pull out corresponding batch covariance
Pvel = P0(4:6, 4:6, n);

% Determine eigen values and vectors
[vecsVel, valsVel] = eig(Pvel);

% Generate the general ellipsoid points for plotting
[ellipVel_x, ellipVel_y, ellipVel_z] = ellipsoid(0, 0, 0, 1, 1, 1);

% Scale the ellipsoid points w.r.t eigenvalues
ellipVel = [ellipVel_x(:), ellipVel_y(:), ellipVel_z(:)]*sqrt(valsVel)*vecsVel';
ellipVel_x(:) = ellipVel(:, 1);
ellipVel_y(:) = ellipVel(:, 2);
ellipVel_z(:) = ellipVel(:, 3);

% Plot the ellipsoid points
%ellipPosFig = figure;
subplot(1, 2, 2)
surf(ellipVel_x, ellipVel_y, ellipVel_z), xlabel('X-Vel Deviation (m/s)'), ...
    ylabel('Y-Vel Deviation (m/s)'), zlabel('Z-Vel Deviation (m/s)'), ...
    title('Velocity Deviation Error Ellipsoid')
hold off


%-------------------------------------------------------------------------%
%                            ECI Orbit Plot                               %
%-------------------------------------------------------------------------%

% Plot the orbit trial
[xE, yE, zE] = sphere;
xE = xE * const.Re; % [m]
yE = yE * const.Re; % [m]
zE = zE * const.Re; % [m]

xOrbit = zeros(1, length(tspan));
yOrbit = zeros(1, length(tspan));
zOrbit = zeros(1, length(tspan));

for m = 1:length(tspan)
    xOrbit(m) = X(1, 1, m);
    yOrbit(m) = X(2, 1, m);
    zOrbit(m) = X(3, 1, m);
end

orbitFig = figure;
plot3(xOrbit, yOrbit, zOrbit), grid, ...
xlabel('\bf I (m)'), ylabel('\bf J (m)'), zlabel('\bf K (m)'), ...
title('ECI Post-Fit Orbit Plot')
hold on
surf(xE, yE, zE, 'EdgeColor', 'none')
alpha 0.2
hold off




end