clc
clear

syms x y z vx vy vz p q u wx wy wz      % States
syms Ft Fx Fy dp dr wtm wtrw tx ty tz T trw tprop m_edf
syms m Jx Jy Jz Jrw ltvc lrw ledf g 

%% Matricies 
% Transformation matrix from body angular velocity to Tait-Bryan rates
W =    [1,      0,          -sin(q)       ;
        0,      cos(p),     cos(q)*sin(p) ;
        0,      -sin(p),    cos(q)*cos(p) ];

Winv = simplify(inv(W));

%rotation about the x-axis 
rx =   [1       0                  0; 
        0      cos(dr)     -sin(dr) ; 
        0      sin(dr)     cos(dr)] ;

%rotation matrix about y-axis 
ry =   [cos(dp)     0       sin(dp) ; 
        0           1             0 ;   
        -sin(dp)    0       cos(dp) ];

%rotation matrix about z-axis (won't be needed in our case)
% rz =   [cos(dy)    -sin(dy)         0; 
%         sin(dy)     cos(dy)         0; 
%         0           0               1];

%Rotation Matrix from Body-to-TVC frames 
Rbt = rx*ry;                      %(rotate about X-axis and then about Y-axis - from Body frame to TVC frame

%Rotation Matrix from TVC-to-Body frames 
Rtb = ry*rx;       %or Rtb = simplify(inv(Rbt))

% Rotation matrix from body to world frame, input: roll, pitch, yaw
R = [cos(q)*cos(u), sin(p)*sin(q)*cos(u)-cos(p)*sin(u), cos(p)*sin(q)*cos(u)+sin(p)*sin(u) ;
     cos(q)*sin(u), sin(p)*sin(q)*sin(u)+cos(p)*cos(u), cos(p)*sin(q)*sin(u)-sin(p)*cos(u) ;
     -sin(q),       sin(p)*cos(q),                      cos(p)*cos(q)                     ];
 

% Matrix of mass moment of inertia
J = [Jx 0  0  ;
     0  Jy 0  ;
     0  0  Jz];


%% Forces

% tx = Ft  * sin(dr); 
% ty = Ft  * sin(dp) * cos(dr); 
tx = Ft * dr; 
ty = Ft * dp;

% Body forces
Tmag = sqrt(tx^2 + ty^2 + Ft^2);
 fb = [tx; ty; Tmag];     %rotate the thrust-vector Ft into body frame 
 
 %fb = [(tx/ltvc) - m_edf*sin(dr); -ty/ltvc - m_edf*sin(dp); Ft];     %rotate the thrust-vector Ft into body frame 
 %Tmag = sqrt(T(1,1)^2 + T(2,1)^2 + T(3,1)^2);


 %% Torques


%tb = [fb(1)*(ltvc+ledf*cos(dr)); fb(2)*(ltvc+ledf*cos(dp)); - trw ];
tb = [fb(1)*(ltvc); fb(2)*(ltvc); trw ];


%% States + Input Vectors 
% State vectors used for derivation
nw = [p q u].';     % Attitude (world frame) 
wb = [wx wy wz].';  % Angular velocity (body frame)
pw = [x y z].';     % Position (world frame)
vb = [vx vy vz].';  % Velocity (body frame)

% Total state vector
X = [nw; wb; pw; vb];                   % full state 
X_att = [nw; wb];                       % attitude only 
X_red = [nw; wb; pw(3); vb(3)];         % Reduced state vector (only attitude and altitude)
X_hor = [ pw(1); pw(2); vb(1); vb(2) ]; % Reduced state vector for horizontal movements
X_roll = [nw(1); wb(1); pw(1); vb(1)];

% Full Input vector (TVC + EDF + Reaction Wheel)
U = [dr; dp; trw; Ft];     
%[pitch_deflection_angle; roll_deflection_angle; edf_angular_speed;reaction_wheel_torque]

% Attitude only input (TVC + EDF) 
U_att = [dr; dp; Ft];


% Input vector for horizontal model
 U_hor = [p; q];

% Roll 
 U_roll = [trw];

 %% Rotational dynamics

nw_dot = Winv * wb; %transform wb to tair-bryon 
wb_dot = inv(J) * (tb - cross(wb, J * wb)  );

 %% Translational dynamics

pw_dot = R * vb;
vb_dot = 1/m * ( fb -  R.' * [0 0 m*g].');

% Translational dynamics in world
vw_dot = 1/m * R * fb - [0 0 m*g].';

%% Combined non-linear model

% Full 12 state vector 
f = [ nw_dot  ;
      wb_dot  ;
      pw_dot  ;
      vb_dot ];

% Reduced non-linear model (8 states)
f_red = [ nw_dot     ;
          wb_dot     ;
          pw_dot(3)  ;
          vb_dot(3) ];
      
% Reduced Attitude only non-linear model (6 states)
f_att = [ nw_dot     ;
          wb_dot     ];
         
      
% Horizontal non-linear model
 f_hor = [ pw_dot(1) ;
           pw_dot(2) ;
           vb_dot(1) ;
           vb_dot(2)];
      
      
 f_roll = [ nw_dot(1) ;
           wb_dot(1)  ;
           pw_dot(1)  ;
           vb_dot(1) ];

 %% Linearization

% Using the Jacobian method, the set of nonlinear system equations are
% linearized around the hover point

 A = jacobian(f_att, X_att);
 B = jacobian(f_att, U);

% Reduced model (only z-axis in position)
A2 = jacobian(f_red, X_red);
B2 = jacobian(f_red, U);

%Reduced model for attitude only 
A5 = jacobian(f_att, X_att);
B5 = jacobian(f_att, U_att);

% Horizontal model (only x- and y-direction)
 A3 = jacobian( f_hor, X_hor );
 B3 = jacobian( f_hor, U_hor );


% Single dimension model (roll/x axis)
 A4 = jacobian( f_roll, X_roll );
 B4 = jacobian( f_roll, U_roll );

 %% Constants + Linearization Point Def 

% The A and B matrixes are now filled with partial derivatives, similar to
% an taylor expansion to approximate a nonlinear function/ODE
% We must insert the state- and input-values at the operating point

% All the states is zero at the hover point
x = 0; y = 0; z = 0; vx = 0; vy = 0; vz = 0; p = 0.0; q = 0; u = 0; wx = 0; wy = 0; wz = 0;

% vehicle Constants
m = 2.300;                      %kg
%Jx = 0.058595;
Jx = 0.027758;
Jy = 0.027758;
%Jy = 0.058595;
Jz = 0.012238;                %Jz = 0.0120276855157049;
Jprop = 0.000174454522116462;
Jrw = 0.00174245619164574;
ltvc = 0.135;                  %COM to TVC (m)
lrw = 0.09525;                  % COM to RW (m)    
ledf = .100;
m_edf = .700;
g= 9.807;


% All input is not zero!
dp=0; dr=0; trw=0; Ft = m*g;

%% System Definition 

% Now the A and B matrixes can be evaluted, yield the full linear model
% around the hover point.
A_sys = vpa(subs(A), 4);
B_sys = vpa(subs(B), 4);
C_sys = eye(6);
D_sys = zeros(6,4);

A_sys = double(A_sys);
B_sys = double(B_sys);
C_sys = double(C_sys);
D_sys = double(D_sys);

% Reduced model 
A_red = vpa(subs(A2), 4);
B_red = vpa(subs(B2), 4);
C_red = eye(8);
D_red = zeros(8,4);
 
 A_red = double(A_red);
 B_red = double(B_red);
 C_red = double(C_red);
 D_red = double(D_red);

%reduced model for attitude only 
A_att = vpa(subs(A5), 4);
B_att = vpa(subs(B5), 4);
C_att = eye(6);
D_att = zeros(6,3);
 
 A_att = double(A_att);
 B_att = double(B_att);
 C_att = double(C_att);
 D_att = double(D_att);

% Horizontal model
 A_hor = double(vpa(subs(A3),4));
 B_hor = double(vpa(subs(B3),4));
 C_hor = double(eye(4));
 D_hor = double(zeros(4,2));

% % Reduced model with integral action states
 G_hov = [ 0 0 0 0 0 0 1 0 ]; % z
%    
 A_int = [A_red; G_hov];
 A_int = [A_int zeros(9,1) ];
 B_int = [B_red; zeros(1,4) ];
 C_int = eye(9);
 D_int = zeros(9,4);

% % Horizontal model with integral action states
 G_pos = [ 1 0 0 0; 
          0 1 0 0 ];
% 
 A_hint = [A_hor; G_pos];
 A_hint = [A_hint zeros(6,2) ];
 B_hint = [B_hor; zeros(2,2) ];
 C_hint = eye(6);
 D_hint = zeros(6,2);

 %% Open Loop dynamics
% 
sys = ss(A_sys,B_sys,C_sys,D_sys);
sys_red = ss(A_red,B_red,C_red,D_red);
sys_att = ss(A_att,B_att,C_att,D_att);
sys_int = ss(A_int,B_int,C_int, D_int);
sys_hor = ss(A_hor, B_hor, C_hor, D_hor);
sys_hint = ss(A_hint, B_hint, C_hint, D_hint);

%% Design controller

%manual values (trial + error)
 qr = 1.2;
 qp = 1.2;
 qy = 0.6;
 qgx = 0.080;
 qgy = 0.080;
 qgz = 3.33;
 qz = .5;
 qvz = 0.06;
 rr = 20;
 rp = 20;
 rrw = 8.71;
 rthrust = 0.01041;   

%Bryson's rule
% qr = 1/0.4^2;
% qp = 1/0.4^2;
% qy = 1/5^2;
% qgx = 1/.8^2;
% qgy = 1/.8^2;
% qgz = 1/3^2;
% qz = 1/2.0^2;
% qvz = 1/4^2;
% rr = 1/0.1^2;
% rp = 1/0.1^2;
% rrw = 1/2^2;
% rthrust = 1/31^2;   



%% Bryson's Rule. 
% Max angle of 0.3 radians. Maximum angular rate of 5 rad/second
Q = [ qr       0        0           0       0       0       0       0       ;  % Roll
      0        qp       0           0       0       0       0       0       ;  % Pitch
      0        0        qy          0       0       0       0       0       ;  % Yaw
      0        0        0           qgx     0       0       0       0       ;  % omega_x
      0        0        0           0       qgy     0       0       0       ;  % omega_y
      0        0        0           0       0       qgz     0       0       ;  % omega_z
      0        0        0           0       0       0       qz      0       ;  % z
      0        0        0           0       0       0       0       qvz      ];  % v_z

Q_red = Q;  
Q_att = Q(1:6,1:6);
Q_sys = Q_att;          %not needed but makes things neater
Q_pos = [.1 0 0 0;
         0 .1 0 0;
         0 0 .25 0;
         0 0 0 0.25];
      

% Max actuation angle of +- 15 degrees 
R = [ rr       0       0        0           ; % dr
      0        rp      0        0           ; % dp
      0        0       rrw        0         ; % wrw
      0        0       0        rthrust    ]; % Ft **
  
R_red = R;              %reduced model R matrix 
R_att = R(1:3,1:3);     %attitude only R matrix 
R_sys = R;          %%not needed but makes things neater
R_pos = [8 0;0 8];
R
Q


%% Optimal Controller 

% Compute "optimal" controller
%K_att = lqr(sys_att, Q_att, R_att);
K_red = lqr(sys_red, Q_red, R_red);

K_sys = lqr(sys, Q_sys, R_sys) ;
K_pos = lqr(sys_hor, Q_pos, R_pos);

%% Closed Loop System Calculation 

% Calcuate closed loop system
%reduced mode (6 states, 4 inputs)
 cl_sys = ss((A_sys - B_sys*K_sys), B_sys, C_sys, D_sys );
 sys_clfb = feedback( sys*K_sys, eye(6));

%attitude only (6 states, 3 inputs)
%cl_sys_att = ss((A_att - B_att*K_att), B_att, C_att, D_att );
%sys_clfb_att7 = feedback( sys_att*K_att, eye(6));

%reduced mode (8 states, 4 inputs)
cl_sys_red = ss((A_red - B_red*K_red), B_red, C_red, D_red );
sys_clfb_red = feedback( sys_red*K_red, eye(8));

K_red
%K_pos = K_pos.*-1;
K_pos
%% 

% Step Response 
%   figure(4)
%   opt = stepDataOptions;                                 %create options object
%   opt.StepAmplitude = 1.26;                                 %set step amplitude 
%   step(sys_clfb_red,  opt, 8)     %run all 3 clfbs on same graph
% 
%  figure(5)
%  X0 = [.08; -.08; 0; 1.5; 1.5; 3.22;2;0];
%  initial(sys_clfb_red, X0, 3);
  