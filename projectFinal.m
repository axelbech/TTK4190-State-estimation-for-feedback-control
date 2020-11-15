% Project in TTK4190 Guidance and Control of Vehicles 
%
% Author:           My name
% Study program:    My study program

clear;
clc;

load('WP.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% USER INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% h  = 0.1;    % sampling time [s]
% Ns = 10000;  % no. of samples

psi_ref = 0 * pi/180;  % desired initial yaw/heading angle (rad)
setPointReached = 0;    % flag used to set new desired yaw/heading
U_d = 7;                % desired cruise speed (m/s)
               
% ship parameters 
m = 17.0677e6;          % mass (kg)
Iz = 2.1732e10;         % yaw moment of inertia about CO (kg m^3)
xg = -3.7;              % CG x-ccordinate (m)
L = 161;                % length (m)
B = 21.8;               % beam (m)
T = 8.9;                % draft (m)
%KT = 0.7;               % propeller coefficient (-)
rho = 1025;             % density of water (kg/m^3)
visc = 1e-6;            % kinematic viscousity at 20 degrees (m/s^2)
eps = 0.001;            % a small number added to ensure that the denominator of Cf is well defined at u=0
k = 0.1;                % form factor giving a viscous correction
t_thr = 0.05;           % thrust deduction number

% propeller parameters
Dia = 3.3;                              % propeller diameter (m)
Ja = 0;                                 % Advance ratio (bollard pull)
PD = 1.5;                               % pitch/diameter ratio
AEAO = 0.65;                            % Blade/area ratio
z = 4;                                  % Number of blades
[KT, KQ] = wageningen(Ja, PD, AEAO, z); % Thrust and torque coefficients

% rudder limitations
delta_max  = 40 * pi/180;        % max rudder angle      (rad)
Ddelta_max = 5  * pi/180;        % max rudder derivative (rad/s)

% added mass matrix about CO
Xudot = -8.9830e5;
Yvdot = -5.1996e6;
Yrdot =  9.3677e5;
Nvdot =  Yrdot;
Nrdot = -2.4283e10;
MA = -[ Xudot 0    0 
        0 Yvdot Yrdot
        0 Nvdot Nrdot ];

% rigid-body mass matrix
MRB = [ m 0    0 
        0 m    m*xg
        0 m*xg Iz ];
    
Minv = inv(MRB + MA); % Added mass is included to give the total inertia

% ocean current in NED
Vc = 1;                             % current speed (m/s)
betaVc = deg2rad(45);               % current direction (rad)

% wind expressed in NED
Vw = 10;                   % wind speed (m/s)
betaVw = deg2rad(135);     % wind direction (rad)
rho_a = 1.247;             % air density at 10 deg celsius
cy = 0.95;                 % wind coefficient in sway
cn = 0.15;                 % wind coefficient in yaw
A_Lw = 10 * L;             % projected lateral area

% linear damping matrix (only valid for zero speed)
T1 = 20; T2 = 20; T6 = 10;

Xu = -(m - Xudot) / T1;
Yv = -(m - Yvdot) / T2;
Nr = -(Iz - Nrdot)/ T6;
D = diag([-Xu -Yv -Nr]);         % zero speed linear damping

% rudder coefficients (Section 9.5)
b = 2;
AR = 8;
CB = 0.8;

lambda = b^2 / AR;
tR = 0.45 - 0.28*CB;
CN = 6.13*lambda / (lambda + 2.25);
aH = 0.75;
xH = -0.4 * L;
xR = -0.5 * L;

X_delta2 = 0.5 * (1 - tR) * rho * AR * CN;
Y_delta = 0.25 * (1 + aH) * rho * AR * CN; 
N_delta = 0.25 * (xR + aH*xH) * rho * AR * CN;   

% input matrix
Bu = @(u_r,delta) [ (1-t_thr)  -u_r^2 * X_delta2 * delta
                        0      -u_r^2 * Y_delta
                        0      -u_r^2 * N_delta            ];
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                    
% Heading Controller
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rudder control law
wb = 0.06;
zeta = 1;
wn = 1 / sqrt( 1 - 2*zeta^2 + sqrt( 4*zeta^4 - 4*zeta^2 + 2) ) * wb;

% linearized sway-yaw model (see (7.15)-(7.19) in Fossen (2021)) used
% for controller design. The code below should be modified.
Yr = 0; % Have we assumed this?
Nv = 0; % Have we assumed this?
N_lin = [-Yv (m-Xudot)*U_d-Yr;
        (Xudot - Yvdot)*U_d-Nv (m*xg-Yrdot)*U_d-Nr];
b_lin = [-2*U_d*Y_delta -2*U_d*N_delta]';
Minv_sway_yaw = Minv(2:3, 2:3);
A_sway_yaw = -Minv_sway_yaw*N_lin;
B_sway_yaw = Minv_sway_yaw*b_lin;
C_sway_yaw = [0 1];
D_sway_yaw = 0;
[num,den] = ss2tf(A_sway_yaw, B_sway_yaw, C_sway_yaw, D_sway_yaw);
% num = 1.0e-04 * [0 0.8638 0.0615]
% den = [1 0.1506 0.0008]
% Laplace domain TF is: r/delta(s) = 1.0e-04*(0.8638s + 0.0615)/(s^2+0.1506s+0.0008)

% 2nd order Nomoto model parameters:
T1_nomoto2 = 181.3575377;
T2_nomoto2 = 6.89246235;
T3_nomoto2 = 14.04552846;
K = 7.6875*10^(-3);
% 1st order Nomoto model parameters:
T_nomoto = T1_nomoto2 + T2_nomoto2 - T3_nomoto2;

% SISO PID Pole placement (algorithm 15.1)
K_p = T_nomoto/K*wn^2;             % Proportional gain
K_d = 2*zeta*wn*T_nomoto/K - 1/K;  % Derivative gain
K_i = (wn/10)*K_p;                 % Integral gain

% initial states
eta = [0 0 180*pi/180]';   % x^n, y^n, yaw/heading (psi) 
nu  = [0.1 0 0]'; % u, v, r (surge, sway, yaw-rate)
nu_dot = [0 0 0]';
delta = 0;
n = 0;
integral_e_psi = 0;
QM = 0;


% 3rd order reference model
wref = 0.1; % reference model natural frequency [rad/s]
xd = [0; 0; 0]; % initial reference model states
Ad = [0 1 0;
      0 0 1;
      -wref^3 -3*wref^2 -3*wref];
Bd = [0 0 wref^3]';

% Initialisation for path following
startWp = WP(:,1);
endWp = WP(:,2);
numWps = length(WP);
wpIdx = 2;
wpDiff = endWp - startWp;
pi_p = atan2(wpDiff(2), wpDiff(1));
% crossTrackErr = crosstrackWpt(endWp(1), endWp(2), startWp(1), startWp(2), eta(1), eta(2));
% dist2Wp = norm(endWp - eta(1:2));
% alongTrackDist = sqrt(dist2Wp^2 - crossTrackErr^2);
% K_p_path = 1 / (alongTrackDist/5 + 500);
yp_int = 0;

% Initialisation for kalman filter
s_r_psi = 0.5 * pi/180; % Standard deviation of heading measurement noise
s_r_r = 0.1 * pi/180; % Standard deviation of heading rate measurement noise
R = [s_r_psi^2 0; 0 s_r_r^2]; % Measurement noise covariance matrix

s_q_psi = s_r_psi; % Standard deviation of heading plant model noise
s_q_r = s_r_r; % Standard deviation of heading rate plant model noise
s_q_b = 0.01; % Standard deviation of rudder bias plant model noise
Q = [s_q_psi 0 0; 0 s_r_r 0; 0 0 s_q_b]; % Plant model noise covariance matrix

Ac = [0 1 0; 0 -1/T_nomoto -K/T_nomoto; 0 0 0];
Bc = [0; K/T_nomoto; 0];
Ec = [0 0; 1 0; 0 1];
Cc = [1; 1; 0];

[~, Bkf] = c2d(Ac,Bc,h);
[Akf, Ekf] = c2d(Ac,Ec,h);
Ckf = Cc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
simdata = zeros(Ns+1,22);                % table of simulation data

for i=1:Ns+1

    t = (i-1) * h;    % time (s)
    
    R = Rzyx(0,0,eta(3));
    
    % current (should be added here)
    u_c = Vc*cos(betaVc - eta(3));
    v_c = Vc*sin(betaVc - eta(3));
    nu_c = [u_c; v_c; 0];
    nu_r = nu - nu_c;
    
    % wind (should be added here)
    u_w = Vw*cos(betaVw - eta(3));
    v_w = Vw*sin(betaVw - eta(3));
    u_rw = nu(1) - u_w;
    v_rw = nu(2) - v_w;
    V_rw = sqrt(u_rw^2 + v_rw^2);
    gamma_rw = -atan2(v_rw, u_rw);
    C_Y = cy*sin(gamma_rw);
    C_N = cn*sin(2*gamma_rw);
    if t > 200
        Ywind = 0.5*rho_a*V_rw^2*C_Y*A_Lw; % expression for wind moment in sway should be added.
        Nwind = 0.5*rho_a*V_rw^2*C_N*A_Lw*L; % expression for wind moment in yaw should be added.
    else
        Ywind = 0;
        Nwind = 0;
    end
    tau_env = [0 Ywind Nwind]';
    
    % state-dependent time-varying matrices
    CRB = m * nu(3) * [ 0 -1 -xg 
                        1  0  0 
                        xg 0  0  ];
                    
    % coriolis due to added mass
    CA = [  0   0   Yvdot * nu_r(2) + Yrdot * nu_r(3)
            0   0   -Xudot * nu_r(1) 
          -Yvdot * nu_r(2) - Yrdot * nu_r(3)    Xudot * nu_r(1)   0];
    N = CRB + CA + D;
    
    % nonlinear surge damping
    Rn = L/visc * abs(nu_r(1));
    Cf = 0.075 / ( (log(Rn) - 2)^2 + eps);
    Xns = -0.5 * rho * (B*L) * (1 + k) * Cf * abs(nu_r(1)) * nu_r(1);
    
    % cross-flow drag
    Ycf = 0;
    Ncf = 0;
    dx = L/10;
    Cd_2D = Hoerner(B,T);
    for xL = -L/2:dx:L/2
        vr = nu_r(2);
        r = nu_r(3);
        Ucf = abs(vr + xL * r) * (vr + xL * r);
        Ycf = Ycf - 0.5 * rho * T * Cd_2D * Ucf * dx;
        Ncf = Ncf - 0.5 * rho * T * Cd_2D * xL * Ucf * dx;
    end
    d = -[Xns Ycf Ncf]';
    
    % Positive 10deg heading setpoint followed by a negative 20deg heading
    % setpoint
%     if eta(3) >= 0.999*psi_ref && setPointReached == 0
%        psi_ref = -180*pi/180;
%        setPointReached = 1;
%     end
%     
    if norm(eta(1:2) - endWp) < L
        if wpIdx ~= numWps
            wpIdx = wpIdx + 1;
            startWp = endWp;
            endWp = WP(:,wpIdx);
            wpDiff = endWp - startWp;
            pi_p = atan2(wpDiff(2), wpDiff(1));
        end
    end     

    [course_ref, yp_int_dot]  = guidance(eta, startWp, endWp, pi_p, yp_int);
    psi_ref = course_ref; % Set heading to desired course
        
    u_d = U_d;
    r_d = 0;
    
    % reference models 
    xd_dot = Ad * xd + Bd * psi_ref;   % Eq. (12.11)
    
    % error signals (psi and r are measurements)
    e_psi = ssa(eta(3) - xd(1));              % yaw angle error (rad)
    e_r = nu(3) - xd(2);                      % yaw rate error (rad/s)
    
    % control law
    delta_c_unsat = -K_p*e_psi - K_i*integral_e_psi - K_d*e_r;              % unsaturated rudder angle command (rad)
  
    % Rudder saturation and dynamics (Sections 9.5.2)
    if abs(delta_c_unsat) >= delta_max
        delta_c = sign(delta_c_unsat)*delta_max;  % Saturation
        integral_e_psi = integral_e_psi - (h/K_i) * (delta_c - delta_c_unsat); % Anti-wind-up
    else
        delta_c = delta_c_unsat;   % No saturation
    end
    
    delta_dot = delta_c - delta;
    if abs(delta_dot) >= Ddelta_max
        delta_dot = sign(delta_dot)*Ddelta_max;  % Angle speed saturation
    end
    
    % propeller and engine dynamics
    Im = 100000; Tm = 10; Km = 0.6;      % propulsion parameters
    Td = (-Xu*(U_d-u_c))/(1-t_thr);
    %Td = -Xu*U_d/(1 - thrustDeductionNumber);   
    %Td = ((m-Xudot)*nu_dot(1)-Xu*(U_d-u_c))/(1-thrustDeductionNumber);
    %Td = -Xu*(U_d-u_c)/(1-thrustDeductionNumber);
    n_d = sign(Td)*sqrt(Td/(rho*Dia^4*KT));             % Desired propeller speed (rps)
    QF = 0;                                                         % Friction torque (assumed zero, not specified)
    Qd = rho * Dia^5 * KQ * abs(n_d) * n_d;                         % Desired propeller moment/torque
    thr = rho * Dia^4 * KT * abs(n) * n;                            % Eq. (9.7) actual thrust
    torque = rho * Dia^5 * KQ * abs(n) * n;                         % Eq. (9.8) actual torque
    Y = Qd/Km;                                                      % Input to engine
    QM_dot = -1/Tm * QM + (Y/Tm)*Km;                                % Main engine dynamics (QM is the produced torque by the engine)
    n_dot = (QM - torque - QF)/Im;                                  % Propeller speed (rps) dynamics
    
    % ship dynamics
    u = [ thr delta ]';
    tau = Bu(nu_r(1),delta) * u;
    nu_dot = [nu(3)*nu_c(2); -nu(3)*nu_c(1); 0] + Minv * (tau_env + tau - N * nu_r - d); % Eq. (10.140)
    eta_dot = R * nu;    % Eq. (10.139)
    integral_e_psi_dot = e_psi; % Augment integral state to system for integral action in controller
    
    beta = asin(nu_r(2) / sqrt( nu_r(1)^2 + nu_r(2)^2 + nu_r(3)^2 )); % Sideslip angle
    beta_c = atan2(nu(2), nu(1)); % Crab angle
    course = eta(3) + beta_c;
    
    psi_noisy = eta(3) + normrnd(0, s_r_psi);
    r_noisy = nu(3) + normrnd(0, s_r_r);
    
    % store simulation data in a table (for testing)
    eta(3) = ssa(eta(3));
    simdata(i,:) = [t n_d delta_c n delta eta' nu' u_d psi_ref r_d beta beta_c, e_psi, integral_e_psi, ssa(course), course_ref, psi_noisy, r_noisy];       
     
    % Euler integration
    eta = euler2(eta_dot,eta,h);
    nu  = euler2(nu_dot,nu,h);
    delta = euler2(delta_dot,delta,h);   
    n  = euler2(n_dot,n,h);    
    xd = euler2(xd_dot,xd,h);
    integral_e_psi = euler2(integral_e_psi_dot,integral_e_psi,h);
    QM = euler2(QM_dot, QM, h);
    yp_int = euler2(yp_int_dot, yp_int, h);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t       = simdata(:,1);                 % s
n_c     = 60 * simdata(:,2);            % rpm
delta_c = (180/pi) * simdata(:,3);      % deg
n       = 60 * simdata(:,4);            % rpm
delta   = (180/pi) * simdata(:,5);      % deg
x       = simdata(:,6);                 % m
y       = simdata(:,7);                 % m
psi     = (180/pi) * simdata(:,8);      % deg
u       = simdata(:,9);                 % m/s
v       = simdata(:,10);                % m/s
r       = (180/pi) * simdata(:,11);     % deg/s
u_d     = simdata(:,12);                % m/s
psi_d   = (180/pi) * simdata(:,13);     % deg
r_d     = (180/pi) * simdata(:,14);     % deg/s
beta    = 180/pi * simdata(:,15);       % deg
beta_c  = 180/pi * simdata(:,16);       % deg
e_psi   = 180/pi * simdata(:,17);
integral_e_psi = 180/pi * simdata(:,18);
course = 180/pi * simdata(:,19);
course_ref = 180/pi * simdata(:,20);
psi_n = 180/pi * simdata(:,21);
r_n = 180/pi * simdata(:,22);

figure(1)
figure(gcf)
subplot(311)
plot(y,x,'linewidth',2); axis('equal')
title('North-East positions (m)');
subplot(312)
plot(t,psi,t,psi_d,'linewidth',2);
title('Actual and desired yaw angles (deg)'); xlabel('time (s)');
subplot(313)
plot(t,r,t,r_d,'linewidth',2);
title('Actual and desired yaw rates (deg/s)'); xlabel('time (s)');

figure(2)
figure(gcf)
subplot(311)
plot(t,u,t,u_d,'linewidth',2);
title('Actual and desired surge velocities (m/s)'); xlabel('time (s)');
subplot(312)
plot(t,n,t,n_c,'linewidth',2);
title('Actual and commanded propeller speed (rpm)'); xlabel('time (s)');
subplot(313)
plot(t,delta,t,delta_c,'linewidth',2);
title('Actual and commanded rudder angles (deg)'); xlabel('time (s)');

figure(3) 
figure(gcf)
subplot(211)
plot(t,u,'linewidth',2);
title('Actual surge velocity (m/s)'); xlabel('time (s)');
subplot(212)
plot(t,v,'linewidth',2);
title('Actual sway velocity (m/s)'); xlabel('time (s)');

figure(4)
figure(gcf)
subplot(211)
plot(t, beta, 'linewidth', 2);
title('Sideslip angle (deg)');
xlabel('time (s)');
subplot(212)
plot(t, beta_c, 'linewidth', 2);
title('Crab angle (deg)');
xlabel('time (s)')

figure(5);
plot(t, course, 'linewidth', 2);
title('Course, desired course and heading (deg)');
xlabel('time (s)')
hold on;
plot(t, course_ref, 'linewidth', 2);
plot(t, psi, 'linewidth', 2);
legend('course', 'desired course', 'heading');
hold off;
% figure(4)
% figure(gcf)
% subplot(211)
% plot(t, e_psi, 'linewidth', 2);
% title('e_{psi} (deg)');
% xlabel('time (s)');
% subplot(212)
% plot(t, integral_e_psi, 'linewidth', 2);
% title('Integral state (deg)');
% xlabel('time (s)')

figure(6)
figure(gcf)
plot(y,x,'linewidth',2); axis('equal')
hold on
siz=size(WP);
for ii=1:(siz(2)-1)   
plot([WP(2,ii), WP(2,ii+1)], [WP(1,ii), WP(1,ii+1)], 'r-x')
end
title('North-East positions (m)');
hold off

figure(7)
subplot(211)
plot(t, psi_n, '--', 'linewidth', 2);
hold on
plot(t, psi, 'linewidth', 2);
hold off
subplot(212)
plot(t, r_n, '--', 'linewidth', 2);
hold on
plot(t, r, 'linewidth', 2);
hold off


