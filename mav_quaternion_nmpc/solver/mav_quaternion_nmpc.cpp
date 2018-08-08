/*
*    This file is part of ACADO Toolkit.
*
*    ACADO Toolkit -- A Toolkit for Automatic Control and Dynamic Optimization.
*    Copyright (C) 2008-2009 by Boris Houska and Hans Joachim Ferreau, K.U.Leuven.
*    Developed within the Optimization in Engineering Center (OPTEC) under
*    supervision of Moritz Diehl. All rights reserved.
*
*    ACADO Toolkit is free software; you can redistribute it and/or
*    modify it under the terms of the GNU Lesser General Public
*    License as published by the Free Software Foundation; either
*    version 3 of the License, or (at your option) any later version.
*
*    ACADO Toolkit is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    Lesser General Public License for more details.
*
*    You should have received a copy of the GNU Lesser General Public
*    License along with ACADO Toolkit; if not, write to the Free Software
*    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*
*/


/**
*    Author David Ariens, Rien Quirynen
*    Date 2009-2013
*    http://www.acadotoolkit.org/matlab 
*/

#include <acado_optimal_control.hpp>
#include <acado_toolkit.hpp>
#include <acado/utils/matlab_acado_utils.hpp>

USING_NAMESPACE_ACADO

#include <mex.h>


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) 
 { 
 
    MatlabConsoleStreamBuf mybuf;
    RedirectStream redirect(std::cout, mybuf);
    clearAllStaticCounters( ); 
 
    mexPrintf("\nACADO Toolkit for Matlab - Developed by David Ariens and Rien Quirynen, 2009-2013 \n"); 
    mexPrintf("Support available at http://www.acadotoolkit.org/matlab \n \n"); 

    if (nrhs != 0){ 
      mexErrMsgTxt("This problem expects 0 right hand side argument(s) since you have defined 0 MexInput(s)");
    } 
 
    TIME autotime;
    DifferentialState p1;
    DifferentialState p2;
    DifferentialState p3;
    DifferentialState v1;
    DifferentialState v2;
    DifferentialState v3;
    DifferentialState e_q1;
    DifferentialState e_q2;
    DifferentialState e_q3;
    DifferentialState e_q4;
    DifferentialState e_w1;
    DifferentialState e_w2;
    DifferentialState e_w3;
    Control f1;
    Control f2;
    Control f3;
    Control f4;
    OnlineData mass; 
    OnlineData inertia1; 
    OnlineData inertia2; 
    OnlineData inertia3; 
    OnlineData drag_coefficient; 
    OnlineData armlength; 
    OnlineData angVel_ref; 
    OnlineData yaw_ref; 
    OnlineData external_forces1; 
    OnlineData external_forces2; 
    OnlineData external_forces3; 
    BMatrix acadodata_M1;
    acadodata_M1.read( "mav_quaternion_nmpc_data_acadodata_M1.txt" );
    BMatrix acadodata_M2;
    acadodata_M2.read( "mav_quaternion_nmpc_data_acadodata_M2.txt" );
    Function acadodata_f1;
    acadodata_f1 << p1;
    acadodata_f1 << p2;
    acadodata_f1 << p3;
    acadodata_f1 << v1;
    acadodata_f1 << v2;
    acadodata_f1 << v3;
    acadodata_f1 << e_q2;
    acadodata_f1 << e_q3;
    acadodata_f1 << e_q4;
    acadodata_f1 << e_w1;
    acadodata_f1 << e_w2;
    acadodata_f1 << e_w3;
    acadodata_f1 << (-1/4.00000000000000000000e+00*9.80659999999999953957e+00*mass+f1);
    acadodata_f1 << (-1/4.00000000000000000000e+00*9.80659999999999953957e+00*mass+f2);
    acadodata_f1 << (-1/4.00000000000000000000e+00*9.80659999999999953957e+00*mass+f3);
    acadodata_f1 << (-1/4.00000000000000000000e+00*9.80659999999999953957e+00*mass+f4);
    Function acadodata_f2;
    acadodata_f2 << p1;
    acadodata_f2 << p2;
    acadodata_f2 << p3;
    acadodata_f2 << v1;
    acadodata_f2 << v2;
    acadodata_f2 << v3;
    OCP ocp1(0, 1.5, 15);
    ocp1.minimizeLSQ(acadodata_M1, acadodata_f1);
    ocp1.minimizeLSQEndTerm(acadodata_M2, acadodata_f2);
    ocp1.subjectTo(5.00000000000000027756e-02 <= f1 <= 6.00000000000000000000e+00);
    ocp1.subjectTo(5.00000000000000027756e-02 <= f2 <= 6.00000000000000000000e+00);
    ocp1.subjectTo(5.00000000000000027756e-02 <= f3 <= 6.00000000000000000000e+00);
    ocp1.subjectTo(5.00000000000000027756e-02 <= f4 <= 6.00000000000000000000e+00);
    DifferentialEquation acadodata_f3;
    acadodata_f3 << dot(p1) == v1;
    acadodata_f3 << dot(p2) == v2;
    acadodata_f3 << dot(p3) == v3;
    acadodata_f3 << dot(v1) == (((((-(f1+f2+f3+f4)*e_q2)*(-e_q4)+(-(f1+f2+f3+f4)*e_q4)*(-e_q2)-(-e_q3)*(f1+f2+f3+f4)*e_q1+(f1+f2+f3+f4)*e_q3*e_q1)*cos(1/2.00000000000000000000e+00*yaw_ref)-((-(f1+f2+f3+f4)*e_q2)*e_q1+(-(f1+f2+f3+f4)*e_q4)*(-e_q3)+(-e_q2)*(f1+f2+f3+f4)*e_q1-(-e_q4)*(f1+f2+f3+f4)*e_q3)*sin(1/2.00000000000000000000e+00*yaw_ref))*cos(1/2.00000000000000000000e+00*yaw_ref)+(((-(f1+f2+f3+f4)*e_q2)*(-e_q4)+(-(f1+f2+f3+f4)*e_q4)*(-e_q2)-(-e_q3)*(f1+f2+f3+f4)*e_q1+(f1+f2+f3+f4)*e_q3*e_q1)*sin(1/2.00000000000000000000e+00*yaw_ref)+((-(f1+f2+f3+f4)*e_q2)*e_q1+(-(f1+f2+f3+f4)*e_q4)*(-e_q3)+(-e_q2)*(f1+f2+f3+f4)*e_q1-(-e_q4)*(f1+f2+f3+f4)*e_q3)*cos(1/2.00000000000000000000e+00*yaw_ref))*(-sin(1/2.00000000000000000000e+00*yaw_ref)))/mass+external_forces1);
    acadodata_f3 << dot(v2) == ((-(((-(f1+f2+f3+f4)*e_q2)*(-e_q4)+(-(f1+f2+f3+f4)*e_q4)*(-e_q2)-(-e_q3)*(f1+f2+f3+f4)*e_q1+(f1+f2+f3+f4)*e_q3*e_q1)*cos(1/2.00000000000000000000e+00*yaw_ref)-((-(f1+f2+f3+f4)*e_q2)*e_q1+(-(f1+f2+f3+f4)*e_q4)*(-e_q3)+(-e_q2)*(f1+f2+f3+f4)*e_q1-(-e_q4)*(f1+f2+f3+f4)*e_q3)*sin(1/2.00000000000000000000e+00*yaw_ref))*(-sin(1/2.00000000000000000000e+00*yaw_ref))+(((-(f1+f2+f3+f4)*e_q2)*(-e_q4)+(-(f1+f2+f3+f4)*e_q4)*(-e_q2)-(-e_q3)*(f1+f2+f3+f4)*e_q1+(f1+f2+f3+f4)*e_q3*e_q1)*sin(1/2.00000000000000000000e+00*yaw_ref)+((-(f1+f2+f3+f4)*e_q2)*e_q1+(-(f1+f2+f3+f4)*e_q4)*(-e_q3)+(-e_q2)*(f1+f2+f3+f4)*e_q1-(-e_q4)*(f1+f2+f3+f4)*e_q3)*cos(1/2.00000000000000000000e+00*yaw_ref))*cos(1/2.00000000000000000000e+00*yaw_ref))/mass+external_forces2);
    acadodata_f3 << dot(v3) == (((-(-(-(f1+f2+f3+f4)*e_q2)*(-e_q2)+(-(f1+f2+f3+f4)*e_q4)*(-e_q4)+(-e_q3)*(f1+f2+f3+f4)*e_q3+(f1+f2+f3+f4)*e_q1*e_q1-9.80659999999999953957e+00*mass)*sin(1/2.00000000000000000000e+00*yaw_ref))*(-sin(1/2.00000000000000000000e+00*yaw_ref))+(-(-(f1+f2+f3+f4)*e_q2)*(-e_q2)+(-(f1+f2+f3+f4)*e_q4)*(-e_q4)+(-e_q3)*(f1+f2+f3+f4)*e_q3+(f1+f2+f3+f4)*e_q1*e_q1-9.80659999999999953957e+00*mass)*cos(1/2.00000000000000000000e+00*yaw_ref)*cos(1/2.00000000000000000000e+00*yaw_ref))/mass+external_forces3);
    acadodata_f3 << dot(e_q1) == (-e_q2*e_w1-e_q3*e_w2-e_q4*e_w3)*5.00000000000000000000e-01;
    acadodata_f3 << dot(e_q2) == (e_q1*e_w1+e_q3*e_w3-e_q4*e_w2)*5.00000000000000000000e-01;
    acadodata_f3 << dot(e_q3) == (e_q1*e_w2-e_q2*e_w3+e_q4*e_w1)*5.00000000000000000000e-01;
    acadodata_f3 << dot(e_q4) == (e_q1*e_w3+e_q2*e_w2-e_q3*e_w1)*5.00000000000000000000e-01;
    acadodata_f3 << dot(e_w1) == ((f2-f4)*armlength-e_w2*e_w3*inertia3+e_w2*inertia2*e_w3)/inertia1;
    acadodata_f3 << dot(e_w2) == (((-f1+f3)*armlength+(angVel_ref+e_w1)*e_w3*inertia3-(angVel_ref+e_w1)*inertia1*e_w3)/inertia2+angVel_ref*e_w3);
    acadodata_f3 << dot(e_w3) == ((-(angVel_ref+e_w1)*e_w2*inertia2+(angVel_ref+e_w1)*inertia1*e_w2+(f1-f2+f3-f4)*drag_coefficient)/inertia3-angVel_ref*e_w2);

    ocp1.setModel( acadodata_f3 );


    ocp1.setNU( 4 );
    ocp1.setNP( 0 );
    ocp1.setNOD( 11 );
    OCPexport ExportModule1( ocp1 );
    ExportModule1.set( GENERATE_MATLAB_INTERFACE, 1 );
    uint options_flag;
    options_flag = ExportModule1.set( HESSIAN_APPROXIMATION, GAUSS_NEWTON );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: HESSIAN_APPROXIMATION");
    options_flag = ExportModule1.set( DISCRETIZATION_TYPE, MULTIPLE_SHOOTING );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: DISCRETIZATION_TYPE");
    options_flag = ExportModule1.set( SPARSE_QP_SOLUTION, FULL_CONDENSING_N2 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: SPARSE_QP_SOLUTION");
    options_flag = ExportModule1.set( INTEGRATOR_TYPE, INT_IRK_GL4 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: INTEGRATOR_TYPE");
    options_flag = ExportModule1.set( NUM_INTEGRATOR_STEPS, 15 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: NUM_INTEGRATOR_STEPS");
    options_flag = ExportModule1.set( QP_SOLVER, QP_QPOASES );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: QP_SOLVER");
    options_flag = ExportModule1.set( HOTSTART_QP, NO );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: HOTSTART_QP");
    options_flag = ExportModule1.set( LEVENBERG_MARQUARDT, 1.000000E-10 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: LEVENBERG_MARQUARDT");
    options_flag = ExportModule1.set( LINEAR_ALGEBRA_SOLVER, GAUSS_LU );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: LINEAR_ALGEBRA_SOLVER");
    options_flag = ExportModule1.set( IMPLICIT_INTEGRATOR_NUM_ITS, 5 );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: IMPLICIT_INTEGRATOR_NUM_ITS");
    options_flag = ExportModule1.set( CG_USE_OPENMP, YES );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: CG_USE_OPENMP");
    options_flag = ExportModule1.set( CG_HARDCODE_CONSTRAINT_VALUES, NO );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: CG_HARDCODE_CONSTRAINT_VALUES");
    options_flag = ExportModule1.set( CG_USE_VARIABLE_WEIGHTING_MATRIX, NO );
    if(options_flag != 0) mexErrMsgTxt("ACADO export failed when setting the following option: CG_USE_VARIABLE_WEIGHTING_MATRIX");
    uint export_flag;
    export_flag = ExportModule1.exportCode( "." );
    if(export_flag != 0) mexErrMsgTxt("ACADO export failed because of the above error(s)!");


    clearAllStaticCounters( ); 
 
} 

