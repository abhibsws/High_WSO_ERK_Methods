#!/usr/bin/env python
# encoding: utf-8
r"""
2D shallow water: radial dam break
==================================

Solve the 2D shallow water equations:

.. math::
    h_t + (hu)_x + (hv)_y = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x + (huv)_y = 0  \\
    (hv)_t + (huv)_x + (hv^2 + \frac{1}{2}gh^2)_y = 0.

The computational domain: [-0.5,0.5]*[-0.5,0.5] and we solve it at t = 0.5.
"""

from __future__ import absolute_import
import numpy as np
from clawpack import riemann
from clawpack.riemann.shallow_roe_with_efix_2D_constants import depth, x_momentum, y_momentum, num_eqn

    
def setup(mx,my,a,b,c,kernel_language='Fortran',use_petsc=False, outdir='./_output',solver_type='sharpclaw',weno_order=5, riemann_solver='roe',disable_output=False,tfinal = 0.5,cfl = 0.8): # added mx,my,a,b,c, solver_type='sharpclaw',weno_order=5,tfinal=1.
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if riemann_solver.lower() == 'roe':
        rs = riemann.shallow_roe_with_efix_2D
    elif riemann_solver.lower() == 'hlle':
        rs = riemann.shallow_hlle_2D

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(rs)
        solver.limiters = pyclaw.limiters.tvd.MC
        solver.dimensional_split=1
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(rs)
        solver.weno_order = weno_order # added 
        solver.time_integrator = 'RK'# added 
        solver.a, solver.b, solver.c = a, b, c # added 
        #solver.dq_src = dq_forcing_2D_SWE # added, this is how the source term is implemented in the sharpclaw solver

    solver.user_bc_lower = time_dep_bc_lower # added
    solver.user_bc_upper = time_dep_bc_upper # added
    solver.bc_lower[0] = pyclaw.BC.custom # left or west boundary
    solver.bc_upper[0] = pyclaw.BC.custom # right or east boundary
    solver.bc_lower[1] = pyclaw.BC.custom # bottom or south boundary
    solver.bc_upper[1] = pyclaw.BC.custom # top or north boundary
    solver.cfl_max = 1.0 # added
    solver.cfl_desired = cfl # added 
        
    # Domain:
    xlower = -0.5; xupper = 0.5
    ylower = -0.5; yupper = 0.5
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])

    state = pyclaw.State(domain,num_eqn)

    # Gravitational constant
    state.problem_data['grav'] = 1.0
    
    # initial condition
    xc = state.grid.x.centers; dx = xc[1]-xc[0]
    yc = state.grid.y.centers; dy = yc[1]-yc[0] 
    X, Y = state.p_centers; t0 = 0 # initial condition

    state.q[depth     ,:,:] = (d_int_h(X+dx/2,Y+dy/2,t0) + d_int_h(X-dx/2,Y-dy/2,t0) - d_int_h(X+dx/2,Y-dy/2,t0) - d_int_h(X-dx/2,Y+dy/2,t0))/(dx*dy)
    state.q[x_momentum,:,:] = (d_int_hu(X+dx/2,Y+dy/2,t0) + d_int_hu(X-dx/2,Y-dy/2,t0) - d_int_hu(X+dx/2,Y-dy/2,t0) - d_int_hu(X-dx/2,Y+dy/2,t0))/(dx*dy)
    state.q[y_momentum,:,:] = (d_int_hv(X+dx/2,Y+dy/2,t0) + d_int_hv(X-dx/2,Y-dy/2,t0) - d_int_hv(X+dx/2,Y-dy/2,t0) - d_int_hv(X-dx/2,Y+dy/2,t0))/(dx*dy)

    claw = pyclaw.Controller()
    claw.tfinal = tfinal
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    if disable_output:
        claw.output_format = None
    claw.outdir = outdir
    claw.num_output_times = 10
    claw.keep_copy = True

    return claw

#---------------------------------------------------------------------------------------#
# Boundary conditions: point values
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
                qbc[0,i,j] = (d_int_h(b1,d1,t) + d_int_h(a1,c1,t) - d_int_h(b1,c1,t) - d_int_h(a1,d1,t))/(dx*dy)
                qbc[1,i,j] = (d_int_hu(b1,d1,t) + d_int_hu(a1,c1,t) - d_int_hu(b1,c1,t) - d_int_hu(a1,d1,t))/(dx*dy)  
                qbc[2,i,j] = (d_int_hv(b1,d1,t) + d_int_hv(a1,c1,t) - d_int_hv(b1,c1,t) - d_int_hv(a1,d1,t))/(dx*dy)
                
    if dim.name == 'y': # bottom edge bc
        for i in range(nx): 
            a2 = (xl-num_ghost*dx)+i*dx; b2 = a2+dx
            for j in range(num_ghost):
                c2 = yl-dy*(num_ghost-j); d2 = c2+dy
                qbc[0,i,j] = (d_int_h(b2,d2,t) + d_int_h(a2,c2,t) - d_int_h(b2,c2,t) - d_int_h(a2,d2,t))/(dx*dy) 
                qbc[1,i,j] = (d_int_hu(b2,d2,t) + d_int_hu(a2,c2,t) - d_int_hu(b2,c2,t) - d_int_hu(a2,d2,t))/(dx*dy)
                qbc[2,i,j] = (d_int_hv(b2,d2,t) + d_int_hv(a2,c2,t) - d_int_hv(b2,c2,t) - d_int_hv(a2,d2,t))/(dx*dy)


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
                qbc[0,nx-num_ghost+i,j] = (d_int_h(b1,d1,t) + d_int_h(a1,c1,t) - d_int_h(b1,c1,t) - d_int_h(a1,d1,t))/(dx*dy)
                qbc[1,nx-num_ghost+i,j] = (d_int_hu(b1,d1,t) + d_int_hu(a1,c1,t) - d_int_hu(b1,c1,t) - d_int_hu(a1,d1,t))/(dx*dy)
                qbc[2,nx-num_ghost+i,j] = (d_int_hv(b1,d1,t) + d_int_hv(a1,c1,t) - d_int_hv(b1,c1,t) - d_int_hv(a1,d1,t))/(dx*dy)  
                       
    if dim.name == 'y': # top edge bc
        for i in range(nx):
            a2 = (xl-num_ghost*dx)+i*dx; b2 = a2+dx
            for j in range(num_ghost): 
                c2 = yu+dy*j; d2 = c2+dy
                qbc[0,i,ny-num_ghost+j] = (d_int_h(b2,d2,t) + d_int_h(a2,c2,t) - d_int_h(b2,c2,t) - d_int_h(a2,d2,t))/(dx*dy) 
                qbc[1,i,ny-num_ghost+j] = (d_int_hu(b2,d2,t) + d_int_hu(a2,c2,t) - d_int_hu(b2,c2,t) - d_int_hu(a2,d2,t))/(dx*dy)
                qbc[2,i,ny-num_ghost+j] = (d_int_hv(b2,d2,t) + d_int_hv(a2,c2,t) - d_int_hv(b2,c2,t) - d_int_hv(a2,d2,t))/(dx*dy)

#---------------------------------------------------------------------------------------#
# Double integral of the manufactured solution and the cell averages of the initial conditions
#---------------------------------------------------------------------------------------#
def d_int_h(x,y,t):
    return x*y*(6.0*t**2 - 0.5*x**2 - 0.5*y**2 + 6.0)/(3.0*t**4 + 6.0*t**2 + 3.0)

def d_int_hu(x,y,t):
    return t*x**2*y*(48.0*t**2 - 6.0*x**2 - 4.0*y**2 + 48.0)/(48.0*t**6 + 144.0*t**4 + 144.0*t**2 + 48.0)

def d_int_hv(x,y,t):
    return t*x*y**2*(24.0*t**2 - 2.0*x**2 - 3.0*y**2 + 24.0)/(24.0*t**6 + 72.0*t**4 + 72.0*t**2 + 24.0)
#---------------------------------------------------------------------------------------#
