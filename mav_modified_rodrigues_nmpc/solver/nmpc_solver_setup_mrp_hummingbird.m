clc;
clear all;
close all;

Ts = 0.1; 
EXPORT = 1;


DifferentialState p(3,1) v(3,1) e_mrp(3,1) e_w(3,1);
Control f1 f2 f3 f4;

OnlineData mass;
OnlineData inertia(3); 
OnlineData drag_coefficient; 
OnlineData armlength;
OnlineData angVel_ref;
OnlineData yaw_ref; 
OnlineData external_forces(3);


n_XD = length(diffStates);
n_U = length(controls);

g = 9.8066; 
m=mass; 
J=diag(inertia); 
inv_J=diag([1/inertia(1), 1/inertia(2), 1/inertia(3)]); 
l= armlength; 
c=drag_coefficient; 

%references
q_ID=[cos(yaw_ref/2); 0; 0; sin(yaw_ref/2)];
u_ss=(m*g/4)*ones(4,1); 
u_ss_total=u_ss(1)+u_ss(2)+u_ss(3)+u_ss(4); 
w_B_ID= [angVel_ref;0;0];  
%% Differential Equation
S=[0 -e_mrp(3) e_mrp(2);e_mrp(3) 0 -e_mrp(1); -e_mrp(2) e_mrp(1) 0];

eq0= ( 1-(e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2) )/(1 + e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 );
eq13= (2/(1 + e_mrp(1)^2 + e_mrp(2)^2 + e_mrp(3)^2 )) * e_mrp;


aux_var=[0;rotate_quat([eq0;eq13],[0;0;0;f1+f2+f3+f4])]-[0; 0; 0; u_ss_total];  
vdot=rotate_quat(q_ID,aux_var);


mq_errdot= 0.25*((1-e_mrp'*e_mrp)*eye(3) + 2*S + 2*e_mrp*e_mrp')*e_w;

M=[(f2-f4)*l; (f3-f1)*l; c*(f1-f2+f3-f4)]; %for hummingbird


f = dot([p; v; e_mrp; e_w]) == ...
    [v;...
    (1/m)* vdot + external_forces; ...
    mq_errdot;...
    inv_J*(M-cross( (e_w+w_B_ID), J*(e_w+w_B_ID) )) + cross(e_w, w_B_ID);...
    ];

h = [p;...
    v;...
    e_mrp;...
    e_w;...
    [f1 f2 f3 f4]' - u_ss];

hN = [p;...
    v];

%% MPCexport
acadoSet('problemname', 'mav_modified_rodrigues_nmpc'); 

N = 15; 
ocp = acado.OCP( 0.0, N*Ts, N );

W_mat = eye(length(h));
WN_mat = eye(length(hN));
W = acado.BMatrix(W_mat);
WN = acado.BMatrix(WN_mat);

ocp.minimizeLSQ( W, h );
ocp.minimizeLSQEndTerm( WN, hN );
ocp.subjectTo(0.05 <= [f1; f2; f3; f4] <= 6.0);
%ocp.subjectTo(-5 <= e_w <= 5);
ocp.setModel(f);


mpc = acado.OCPexport( ocp );
mpc.set( 'HESSIAN_APPROXIMATION',       'GAUSS_NEWTON'      );
mpc.set( 'DISCRETIZATION_TYPE',         'MULTIPLE_SHOOTING' );
mpc.set( 'SPARSE_QP_SOLUTION',        'FULL_CONDENSING_N2'  ); 
mpc.set( 'INTEGRATOR_TYPE',             'INT_IRK_GL4'       );
mpc.set( 'NUM_INTEGRATOR_STEPS',         N                  );
mpc.set( 'QP_SOLVER',                   'QP_QPOASES'    	);
mpc.set( 'HOTSTART_QP',                 'NO'             	);
mpc.set( 'LEVENBERG_MARQUARDT',          1e-10				);


mpc.set( 'LINEAR_ALGEBRA_SOLVER',        'GAUSS_LU'         ); 
mpc.set( 'IMPLICIT_INTEGRATOR_NUM_ITS',  5                  );
mpc.set( 'CG_USE_OPENMP',                'YES'              );
mpc.set( 'CG_HARDCODE_CONSTRAINT_VALUES', 'NO'              );
mpc.set( 'CG_USE_VARIABLE_WEIGHTING_MATRIX', 'NO'           );


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

