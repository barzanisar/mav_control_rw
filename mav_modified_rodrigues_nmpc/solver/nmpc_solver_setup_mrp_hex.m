clc;
clear all;
close all;

Ts = 0.1; %prediction sampling time %0.05
EXPORT = 1;

DifferentialState e_p(3,1) e_v(3,1) e_mrp(3,1) e_w(3,1);
Control f1 f2 f3 f4 f5 f6;

OnlineData mass;
OnlineData inertia(3); % vector Ixx, Iyy, Izz 1x3 vec?
OnlineData drag_coefficient; % c=0.0001;
OnlineData armlength;
OnlineData angVel_ref;
OnlineData yaw_ref; %pi/4  //should reference and arm length be online data?
OnlineData external_forces(3);
%OnlineData quat_norm_gain;

n_XD = length(diffStates);
n_U = length(controls);

g = 9.8066; %g = 9.81;
m=mass; %1.9;
J=diag(inertia); %diag([0.09,0.09,0.17]); %Inertia
inv_J=diag([1/inertia(1), 1/inertia(2), 1/inertia(3)]); %inv(J);
l= armlength; %length
c=drag_coefficient; %0.0001; %drag coefficient

% %references
% yaw_ref=pi/4;
q_ID=[cos(yaw_ref/2); 0; 0; sin(yaw_ref/2)];
%R_ID=[cos(yaw_ref) -sin(yaw_ref) 0; sin(yaw_ref) cos(yaw_ref) 0; 0 0 1];
u_ss=(m*g/n_U)*ones(n_U,1); %from outside 
u_ss_total= m*g; %mg
w_B_ID= [angVel_ref;0;0];  % should be online data since desired angular vel??
%% Differential Equation

S=[0 -e_mrp(3) e_mrp(2);e_mrp(3) 0 -e_mrp(1); -e_mrp(2) e_mrp(1) 0];

eq0= ( 1-(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2) )/(1 + e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 );
eq13= (2/(1 + e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 )) * e_mrp;

A=[0.0435778713738291,-0.0871557427476582,0.0435778713738291,0.0435778713738291,-0.0871557427476582,0.0435778713738291;-0.0754790873051733,0,0.0754790873051733,-0.0754790873051733,0,0.0754790873051733;0.996194698091746,0.996194698091746,0.996194698091746,0.996194698091746,0.996194698091746,0.996194698091746;0.143469078891783,0.286938157783566,0.143469078891783,-0.143469078891783,-0.286938157783566,-0.143469078891783;-0.248495733955676,0,0.248495733955676,0.248495733955676,0,-0.248495733955676;-0.0419218334972761,0.0419218334972761,-0.0419218334972761,0.0419218334972761,-0.0419218334972761,0.0419218334972761];
F_M= A*[f1;f2;f3;f4;f5;f6];

%err_thrust=[0;rotate_quat([eq0;eq13],[0;F_M(1:3)])]-[0; 0; 0; u_ss_total];   %-[0; 0; 0; u_ss_total];
err_thrust=[0;rotate_quat([eq0;eq13],[0;0;0;F_M(3)])]-[0; 0; 0; u_ss_total];   %-[0; 0; 0; u_ss_total];
vdot=rotate_quat(q_ID,err_thrust);

% in=[            (e_mrp(1)^2 + 1)/(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 + 1), -(e_mrp(3) - e_mrp(1)*e_mrp(2))/(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 + 1),  (e_mrp(2) + e_mrp(1)*e_mrp(3))/(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 + 1); ...
%   (e_mrp(3) + e_mrp(1)*e_mrp(2))/(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 + 1),            (e_mrp(2)^2 + 1)/(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 + 1), -(e_mrp(1) - e_mrp(2)*e_mrp(3))/(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 + 1);...
%  -(e_mrp(2) - e_mrp(1)*e_mrp(3))/(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 + 1),  (e_mrp(1) + e_mrp(2)*e_mrp(3))/(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 + 1),            (e_mrp(3)^2 + 1)/(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 + 1)];
%   
%  
% R_DB= (eye(3) + S)*(eye(3) + S)*in * in ;
% err_thrust=R_DB * [0;0;f1+f2+f3+f4]-[0; 0; u_ss_total];
%vdot=R_ID*err_thrust;


mq_errdot= 0.25*((1-e_mrp'*e_mrp)*eye(3) + 2*S + 2*e_mrp*e_mrp')*e_w;

M = F_M(4:6); %moments % include sign of ftotal in yawing moment
%M=[0;0;F_M(6)];


f = dot([e_p; e_v; e_mrp; e_w]) == ...
    [e_v;...
    (1/m)* vdot + external_forces; ...
    mq_errdot;... %+0.5*r*e_q*(1/(e_q'*e_q)-1) %baumgarte stabilisation
    inv_J*(M-cross( (e_w+w_B_ID), J*(e_w+w_B_ID) )) + cross(e_w, w_B_ID);...
    ];

h = [e_p;...
    e_v;...
    e_mrp;...
    e_w;...
    [f1 f2 f3 f4 f5 f6]' - u_ss];

hN = [e_p;...
    e_v];

%% MPCexport
acadoSet('problemname', 'mav_modified_rodrigues_nmpc'); %'barza_mpc'

N = 20; %40
ocp = acado.OCP( 0.0, N*Ts, N );

W_mat = eye(length(h));
WN_mat = eye(length(hN));
W = acado.BMatrix(W_mat);
WN = acado.BMatrix(WN_mat);

ocp.minimizeLSQ( W, h );
ocp.minimizeLSQEndTerm( WN, hN );
ocp.subjectTo(-8.0 <= [f1; f2; f3; f4; f5; f6] <= 8.0);
%ocp.subjectTo(-5 <= e_w <= 5);
ocp.subjectTo(-2.5 <= e_mrp <= 2.5);
%ocp.subjectTo(cos(65*pi/180) <= R_DB(3,3)  );
ocp.setModel(f);


mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',        'FULL_CONDENSING_N2'  ); %FULL_CONDENsinG_N2
mpc.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL4'       );
mpc.set( 'NUM_INTEGRATOR_STEPS',         N                  );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);
mpc.set( 'HOTSTART_QP',                 'NO'             	);
mpc.set( 'LEVENBERG_MARQUARDT',          1e-10				);


mpc.set( 'LINEAR_ALGEBRA_SOLVER',        'GAUSS_LU'         ); %Do we need this?
mpc.set( 'IMPLICIT_INTEGRATOR_NUM_ITS',  5                  );
mpc.set( 'CG_USE_OPENMP',                'YES'              );
mpc.set( 'CG_HARDCODE_CONSTRAINT_VALUES', 'NO'              );
mpc.set( 'CG_USE_VARIABLE_WEIGHTING_MATRIX', 'NO'           );

% if EXPORT
%     mpc.exportCode( 'export_MPC' );
%     copyfile('../../../../../../external_packages/qpoases', 'export_MPC/qpoases', 'f')
%     
%     cd export_MPC
%     make_acado_solver('../acado_MPCstep')
%     cd ..
% end

if EXPORT
    mpc.exportCode('.');
end


function [rotated_quat]=rotate_quat(q,v) 
% q and v are 4x1 quats
anss= quat_mult(quat_mult(q,v), [q(1); -q(2:4)]);
rotated_quat=anss(2:4); %to covert to 3x1 vec
end

function [mult_quat]=quat_mult(q,p)
%q and p are 4x1 quats
mult_quat=[ p(1)*q(1) - p(2)*q(2) - p(3)*q(3) - p(4)*q(4), p(1)*q(2) + p(2)*q(1) - p(3)*q(4) + p(4)*q(3), p(1)*q(3) + p(3)*q(1) + p(2)*q(4) - p(4)*q(2), p(1)*q(4) - p(2)*q(3) + p(3)*q(2) + p(4)*q(1)]';
%returns 4x1 quat 
end

