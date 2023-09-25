#!/usr/bin/env python
# encoding: utf-8

r"""
Shallow water flow with manufactured solution.
=============================================

Solve the one-dimensional shallow water equations:

.. math::
    h_t + (hu)_x & = f_u \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x & = f_hu.

on the domain [0,1]*(0,0.5]. Here h is the depth, u is the velocity, and g is the gravitational constant.
"""

from __future__ import absolute_import
import numpy as np
from clawpack import riemann
from clawpack.riemann.shallow_roe_with_efix_1D_constants import depth, momentum, num_eqn
import scipy.integrate as integrate

def setup(mx, a, b, c, use_petsc=False,kernel_language='Fortran',outdir='./_output',solver_type='sharpclaw', weno_order=5, 
           riemann_solver='roe', disable_output=False, tfinal=0.5, cfl = 0.8): 
    # Inputs: mx,a,b,c.
    # solver_type='sharpclaw', weno_order=5, tfinla = 0.5, cfl = 0.8. These are default parameters.

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw # 

    if kernel_language == 'Python':
        if riemann_solver.lower() == 'roe':
            raise Exception('Python Roe solver not implemented.')
        elif riemann_solver.lower() == 'hlle':
            rs = riemann.shallow_1D_py.shallow_hll_1D
    elif kernel_language == 'Fortran': #
        if riemann_solver.lower() == 'roe':
            rs = riemann.shallow_roe_with_efix_1D
        elif riemann_solver.lower() == 'hlle':
            rs = riemann.shallow_hlle_1D
 
    solver = pyclaw.SharpClawSolver1D(rs)
    solver.weno_order = weno_order # added 
    solver.time_integrator = 'RK'# added 
    solver.a, solver.b, solver.c = a, b, c # added 
    solver.dq_src = dq_forcing_SWE # added. The source term is implemented in the sharpclaw solver
  

    solver.kernel_language = kernel_language
    
    solver.user_bc_lower = time_dep_bc_lower # added
    solver.user_bc_upper = time_dep_bc_upper # added
    solver.bc_lower[0] = pyclaw.BC.custom # added: left boundary
    solver.bc_upper[0] = pyclaw.BC.custom # added: right boundary

    solver.cfl_max = 2.0 # added 
    solver.cfl_desired = cfl # added 

    xlower = 0
    xupper = 1
    # mx = 500 # added as input
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,num_eqn)

    # Gravitational constant
    state.problem_data['grav'] = 1.0
    state.problem_data['dry_tolerance'] = 1e-3
    state.problem_data['sea_level'] = 0.0
    
    # Initial data # added
    xc = state.grid.x.centers; dx = xc[1]-xc[0]; t0 = 0
    for i in range(len(xc)):
        I1 = integrate.quad(lambda x: fun_h(x,t0), xc[i]-dx/2,xc[i]+dx/2)
        I2 = integrate.quad(lambda x: fun_hu(x,t0), xc[i]-dx/2,xc[i]+dx/2)
        state.q[depth,i] = I1[0]/dx
        state.q[momentum,i] = I2[0]/dx


    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.tfinal = tfinal # 
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    #claw.setplot = setplot

    return claw

#-----------------------------------------------------------------------#
# Required functions: to change the manufactured solution only change fun_h and fun_u
def fun_h(x,t):
    return (1+x)/(1+t) 

def fun_u(x,t):
    return (1+x**2)/(0.5+t) 

# momentum function
def fun_hu(x,t):
    return fun_h(x,t)*fun_u(x,t)

# Forcing functions for mass and momentum: change these functions depending on the manufactured solution from the other
# file: Symbolic_Cal_Forcing_SWE.py
def fun_forcing_h(x,t):
    return 2*x*(x + 1)/((t + 0.5)*(t + 1)) - (x + 1)/(t + 1)**2 + (x**2 + 1)/((t + 0.5)*(t + 1)) 

def fun_forcing_hu(x,t):
    return 4*x*(x + 1)*(x**2 + 1)/((t + 0.5)**2*(t + 1)) + 0.5*(2*x + 2)/(t + 1)**2 - (x + 1)*(x**2 + 1)/((t + 0.5)*(t + 1)**2) - (x + 1)*(x**2 + 1)/((t + 0.5)**2*(t + 1)) + (x**2 + 1)**2/((t + 0.5)**2*(t + 1)) 


#------------------------------------------------------------------------------------------------#
def int_fun_forcing_h(x,t):
    return (3.0*x**3/(3.0*t**2 + 4.5*t + 1.5) + x**2*(1.0*t + 1.5)/(2.0*t**3 + 5.0*t**2 + 4.0*t + 1.0) + 0.5*x/(1.0*t**3 + 2.5*t**2 + 2.0*t + 0.5))

def int_fun_forcing_hu(x,t):
    return (5.0*x**5/(5.0*t**3 + 10.0*t**2 + 6.25*t + 1.25) + x**4*(2.0*t + 2.5)/(4.0*t**4 + 12.0*t**3 + 13.0*t**2 + 6.0*t + 1.0) + x**3*(4.0*t + 4.5)/(3.0*t**4 + 9.0*t**3 + 9.75*t**2 + 4.5*t + 0.75) + x**2*(1.0*t**2 + 3.0*t + 2.75)/(2.0*t**4 + 6.0*t**3 + 6.5*t**2 + 3.0*t + 0.5) + x*(1.0*t - 0.5)/(1.0*t**3 + 2.5*t**2 + 2.0*t + 0.5))

#------------------------------------------------------------------------------------------------#

# These functions take cell averages which are required for higher order convergence: Cell averages
#------------------------------------------------------------------------------------------------#
def dq_forcing_SWE(solver,state,dt):
    xc = state.grid.x.centers; dx = xc[1]-xc[0]; t=state.t
    dq = np.empty(state.q.shape)
    for i in range(len(xc)):
        dq[0,i] = dt*(int_fun_forcing_h(xc[i]+dx/2,t)-int_fun_forcing_h(xc[i]-dx/2,t))/dx
        dq[1,i] = dt*(int_fun_forcing_hu(xc[i]+dx/2,t)-int_fun_forcing_hu(xc[i]-dx/2,t))/dx  
    return dq

def time_dep_bc_lower(state,dim,t,qbc,auxbc,num_ghost):
    xl = state.grid.x.lower; xu = state.grid.x.upper; xc = state.grid.x.centers; dx = xc[1]-xc[0]
    nx = len(xc)+2*num_ghost
    for i in range(num_ghost):
        a1 = xl-dx*(num_ghost-i); b1 = xl-dx*(num_ghost-(i+1))
        I1 = integrate.quad(lambda x: fun_h(x,t),a1,b1)
        I2 = integrate.quad(lambda x: fun_hu(x,t),a1,b1)
        qbc[0,i] = I1[0]/dx
        qbc[1,i] = I2[0]/dx
                    
def time_dep_bc_upper(state,dim,t,qbc,auxbc,num_ghost):
    xl = state.grid.x.lower; xu = state.grid.x.upper; xc = state.grid.x.centers; dx = xc[1]-xc[0]
    nx = len(xc)+2*num_ghost
    for i in range(num_ghost):
        a1 = xu+dx*i; b1 = xu+dx*(i+1);
        I1 = integrate.quad(lambda x: fun_h(x,t),a1,b1)
        I2 = integrate.quad(lambda x: fun_hu(x,t),a1,b1)
        qbc[0,num_ghost+len(xc)+i] = I1[0]/dx
        qbc[1,num_ghost+len(xc)+i] = I2[0]/dx              
#------------------------------------------------------------------------------------------------#  