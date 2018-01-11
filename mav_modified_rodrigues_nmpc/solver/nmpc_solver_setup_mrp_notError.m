clc;
clear all;
close all;

Ts = 0.1; 
EXPORT = 1;


DifferentialState p(3,1) v(3,1) mrp(3,1) w(3,1);
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

% %references
q_ID=[cos(yaw_ref/2); 0; 0; sin(yaw_ref/2)];
%R_ID=[cos(yaw_ref) -sin(yaw_ref) 0; sin(yaw_ref) cos(yaw_ref) 0; 0 0 1];
u_ss=(m*g/4)*ones(4,1); 
u_ss_total=m*g; 
%w_B_ID= [angVel_ref;0;0];  
%% Differential Equation
S=[0 -mrp(3) mrp(2);mrp(3) 0 -mrp(1); -mrp(2) mrp(1) 0];

q0= ( 1-(mrp(1)^2 + mrp(2)^2 + mrp(3)^2) )/(1 + mrp(1)^2 + mrp(2)^2 + mrp(3)^2 );
q13= (2/(1 + mrp(1)^2 + mrp(2)^2 + mrp(3)^2 )) * mrp;


aux=[0;rotate_quat([q0;q13],[0;0;0;f1+f2+f3+f4])]-[0; 0; 0; u_ss_total];   %-[0; 0; 0; u_ss_total];
vdot=aux(2:4);

e_q=quat_mult([q_ID(1); -q_ID(2:4)],[q0; q13]);
e_mrp= e_q(2:4)/(1+e_q(1));

mrp_dot= 0.25*((1-mrp'*mrp)*eye(3) + 2*S + 2*mrp*mrp')*w;

M=[(f2-f4)*l; (f3-f1)*l; c*(f1-f2+f3-f4)]; 


f = dot([p; v; mrp; w]) == ...
    [v;...
    (1/m)* vdot + external_forces; ... 
    mrp_dot;... 
    inv_J*(M-cross( w, J*w )) ;...
    ];

h = [p;...
    v;...
    e_mrp;...
    w;...
    [f1 f2 f3 f4]' - u_ss]; %- u_ss

hN = [p;...
    v];
%% MPCexport
acadoSet('problemname', 'mav_modified_rodrigues_nmpc'); 

N = 15; %40
ocp = acado.OCP( 0.0, N*Ts, N );

W_mat = eye(length(h));
WN_mat = eye(length(hN));
W = acado.BMatrix(W_mat);
WN = acado.BMatrix(WN_mat);

ocp.minimizeLSQ( W, h );
ocp.minimizeLSQEndTerm( WN, hN );
ocp.subjectTo(0.05 <= [f1; f2; f3; f4] <= 6.0);
%ocp.subjectTo(-5 <= e_w <= 5);
%ocp.subjectTo(cos(65*pi/180) <= R_DB(3,3)  );
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

