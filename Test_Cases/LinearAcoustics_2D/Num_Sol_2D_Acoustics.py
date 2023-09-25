#!/usr/bin/env python
# encoding: utf-8
r"""
Two-dimensional acoustics
=========================

Solve the (linear) acoustics equations with forcing:

.. math:: 
    p_t + K (u_x + v_y) & = f^p \\ 
    u_t + p_x / \rho & = f^u \\
    v_t + p_y / \rho & = f^v.

Here p is the pressure, (u,v) is the velocity, K is the bulk modulus,
and :math:`\rho` is the density. We solve the problem on domain 
[0,1]*[0,1] at a final time T = 1. \rho = 1 and K = 4.
"""
from __future__ import absolute_import
from clawpack import riemann
import numpy as np

def setup(mx, my, a, b, c, kernel_language='Fortran', use_petsc=False, outdir='./_output', 
              solver_type='sharpclaw', weno_order=5, ptwise=False,
              disable_output=False, tfinal = 0.1, cfl = 0.45):
    """
    Example python script for solving the 2d acoustics equations.
    """
    if use_petsc:
        from clawpack import petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver=pyclaw.SharpClawSolver2D(riemann.acoustics_2D)
        solver.weno_order = weno_order # added 
        solver.time_integrator = 'RK'# added 
        solver.a, solver.b, solver.c = a, b, c # added 
        solver.dq_src = dq_forcing_2D_Acoustics # added, this is how the source term is implemented in the sharpclaw solver
       
    solver.user_bc_lower = time_dep_bc_lower # added
    solver.user_bc_upper = time_dep_bc_upper # added
    solver.bc_lower[0] = pyclaw.BC.custom # left or west boundary
    solver.bc_upper[0] = pyclaw.BC.custom # right or east boundary
    solver.bc_lower[1] = pyclaw.BC.custom # bottom or south boundary
    solver.bc_upper[1] = pyclaw.BC.custom # top or north boundary
    solver.cfl_max = 0.5 # added
    solver.cfl_desired = cfl # added 
    
    # Domain
    x = pyclaw.Dimension(0,1.0,mx,name='x')
    y = pyclaw.Dimension(0,1.0,my,name='y')
    domain = pyclaw.Domain([x,y])

    num_eqn = 3
    state = pyclaw.State(domain,num_eqn)

    rho  = 1.0  # Material density
    bulk = 4.0  # Material bulk modulus
    cc = np.sqrt(bulk/rho)  # sound speed
    zz = rho*cc             # impedance
    state.problem_data['rho']= rho
    state.problem_data['bulk']=bulk
    state.problem_data['zz']= zz
    state.problem_data['cc']=cc

    solver.dt_initial=np.min(domain.grid.delta)/state.problem_data['cc']*solver.cfl_desired

    qinit(state)

    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.num_output_times = 10
    claw.tfinal = tfinal


    return claw


#---------------------------------------------------------------------------------------#
# Cell averages of the initial conditions
#---------------------------------------------------------------------------------------#
def qinit(state):
    xc = state.grid.x.centers; dx = xc[1]-xc[0]
    yc = state.grid.y.centers; dy = yc[1]-yc[0] 
    X, Y = state.p_centers; t0 = 0 # initial condition

    state.q[0,:,:] = (d_int_p(X+dx/2,Y+dy/2,t0) + d_int_p(X-dx/2,Y-dy/2,t0) - d_int_p(X+dx/2,Y-dy/2,t0) - d_int_p(X-dx/2,Y+dy/2,t0))/(dx*dy)
    state.q[1,:,:] = (d_int_u(X+dx/2,Y+dy/2,t0) + d_int_u(X-dx/2,Y-dy/2,t0) - d_int_u(X+dx/2,Y-dy/2,t0) - d_int_u(X-dx/2,Y+dy/2,t0))/(dx*dy)
    state.q[2,:,:] = (d_int_v(X+dx/2,Y+dy/2,t0) + d_int_v(X-dx/2,Y-dy/2,t0) - d_int_v(X+dx/2,Y-dy/2,t0) - d_int_v(X-dx/2,Y+dy/2,t0))/(dx*dy)
#---------------------------------------------------------------------------------------#
# Dirichlet boundary conditions: point values
#---------------------------------------------------------------------------------------#
def time_dep_bc_lower(state,dim,t,qbc,auxbc,num_ghost):
    xl = state.grid.x.lower; xu = state.grid.x.upper
    yl = state.grid.y.lower; yu = state.grid.y.upper
    xc = state.grid.x.centers; dx = xc[1]-xc[0]; yc = state.grid.y.centers; dy = yc[1]-yc[0]; 
    nx = len(xc)+2*num_ghost; ny = len(yc)+2*num_ghost
    if dim.name == 'x': # left edge bc
        for i in range(num_ghost):
            a1 = xl-dx*(num_ghost-i); b1 = a1+dx; 
            for j in range(ny): 
                c1 = (yl-num_ghost*dy)+j*dy; d1 = c1+dy
                qbc[0,i,j] = (d_int_p(b1,d1,t) + d_int_p(a1,c1,t) - d_int_p(b1,c1,t) - d_int_p(a1,d1,t))/(dx*dy)
                qbc[1,i,j] = (d_int_u(b1,d1,t) + d_int_u(a1,c1,t) - d_int_u(b1,c1,t) - d_int_u(a1,d1,t))/(dx*dy)  
                qbc[2,i,j] = (d_int_v(b1,d1,t) + d_int_v(a1,c1,t) - d_int_v(b1,c1,t) - d_int_v(a1,d1,t))/(dx*dy)
                
    if dim.name == 'y': # bottom edge bc
        for i in range(nx): 
            a2 = (xl-num_ghost*dx)+i*dx; b2 = a2+dx
            for j in range(num_ghost):
                c2 = yl-dy*(num_ghost-j); d2 = c2+dy
                qbc[0,i,j] = (d_int_p(b2,d2,t) + d_int_p(a2,c2,t) - d_int_p(b2,c2,t) - d_int_p(a2,d2,t))/(dx*dy) 
                qbc[1,i,j] = (d_int_u(b2,d2,t) + d_int_u(a2,c2,t) - d_int_u(b2,c2,t) - d_int_u(a2,d2,t))/(dx*dy)
                qbc[2,i,j] = (d_int_v(b2,d2,t) + d_int_v(a2,c2,t) - d_int_v(b2,c2,t) - d_int_v(a2,d2,t))/(dx*dy)


def time_dep_bc_upper(state,dim,t,qbc,auxbc,num_ghost):
    xl = state.grid.x.lower; xu = state.grid.x.upper
    yl = state.grid.y.lower; yu = state.grid.y.upper
    xc = state.grid.x.centers; dx = xc[1]-xc[0]; yc = state.grid.y.centers; dy = yc[1]-yc[0]; 
    ny = len(yc)+2*num_ghost; nx = len(xc)+2*num_ghost
    if dim.name == 'x':  # right edge bc           
        for i in range(num_ghost):
            a1 = xu+dx*i; b1 = a1+dx
            for j in range(ny): 
                c1 = (yl-num_ghost*dy)+j*dy; d1 = c1+dy
                qbc[0,nx-num_ghost+i,j] = (d_int_p(b1,d1,t) + d_int_p(a1,c1,t) - d_int_p(b1,c1,t) - d_int_p(a1,d1,t))/(dx*dy)
                qbc[1,nx-num_ghost+i,j] = (d_int_u(b1,d1,t) + d_int_u(a1,c1,t) - d_int_u(b1,c1,t) - d_int_u(a1,d1,t))/(dx*dy)
                qbc[2,nx-num_ghost+i,j] = (d_int_v(b1,d1,t) + d_int_v(a1,c1,t) - d_int_v(b1,c1,t) - d_int_v(a1,d1,t))/(dx*dy) 
                   
                    
    if dim.name == 'y': # top edge bc
        for i in range(nx):
            a2 = (xl-num_ghost*dx)+i*dx; b2 = a2+dx
            for j in range(num_ghost): 
                c2 = yu+dy*j; d2 = c2+dy
                qbc[0,i,ny-num_ghost+j] = (d_int_p(b2,d2,t) + d_int_p(a2,c2,t) - d_int_p(b2,c2,t) - d_int_p(a2,d2,t))/(dx*dy) 
                qbc[1,i,ny-num_ghost+j] = (d_int_u(b2,d2,t) + d_int_u(a2,c2,t) - d_int_u(b2,c2,t) - d_int_u(a2,d2,t))/(dx*dy)
                qbc[2,i,ny-num_ghost+j] = (d_int_v(b2,d2,t) + d_int_v(a2,c2,t) - d_int_v(b2,c2,t) - d_int_v(a2,d2,t))/(dx*dy)
#---------------------------------------------------------------------------------------#



#---------------------------------------------------------------------------------------#    
# Change here for different manufactured solutions
#---------------------------------------------------------------------------------------#
# Double integral of the manufactured solution and the cell averages of solutions
#---------------------------------------------------------------------------------------#
# Manufactures solution 
#-------------------------#
def d_int_p(x,y,t):
    return x*y**2/(2*t + 2) + y*(x**2 + 2*x)/(2*t + 2)

def d_int_u(x,y,t):
    return t*x**2*y/(2*t + 2)

def d_int_v(x,y,t):
    return t*x*y**2/(2*t + 2)
# Change here for different manufactured solutions
#---------------------------------------------------------------------------------------#
# Double integral of the forcing functions
#---------------------------------------------------------------------------------------#
def d_int_f_p(x,y,t):
    return -x*y**2/(2*t**2 + 4*t + 2) + y*(16*t**2*x + 16*t*x - x**2 - 2*x)/(2*t**2 + 4*t + 2)

def d_int_f_u(x,y,t):
    return y*(x**2/(2*t**2 + 4*t + 2) + x/(t + 1))

def d_int_f_v(x,y,t):
    return x*y**2/(2*t**2 + 4*t + 2) + x*y/(t + 1)
#---------------------------------------------------------------------------------------#
# Source term in 2D: cell averages of the forcing functions
#---------------------------------------------------------------------------------------#
def dq_forcing_2D_Acoustics(solver,state,dt):
    xc = state.grid.x.centers; dx = xc[1]-xc[0]; t=state.t
    yc = state.grid.y.centers; dy = yc[1]-yc[0];  
    X, Y = state.p_centers;
    dq = np.empty(state.q.shape)
    dq[0,:,:] = dt*(d_int_f_p(X+dx/2,Y+dy/2,t) + d_int_f_p(X-dx/2,Y-dy/2,t) - d_int_f_p(X+dx/2,Y-dy/2,t) - d_int_f_p(X-dx/2,Y+dy/2,t))/(dx*dy)
    dq[1,:,:] = dt*(d_int_f_u(X+dx/2,Y+dy/2,t) + d_int_f_u(X-dx/2,Y-dy/2,t) - d_int_f_u(X+dx/2,Y-dy/2,t) - d_int_f_u(X-dx/2,Y+dy/2,t))/(dx*dy)
    dq[2,:,:] = dt*(d_int_f_v(X+dx/2,Y+dy/2,t) + d_int_f_v(X-dx/2,Y-dy/2,t) - d_int_f_v(X+dx/2,Y-dy/2,t) - d_int_f_v(X-dx/2,Y+dy/2,t))/(dx*dy)

    return dq
#---------------------------------------------------------------------------------------#


