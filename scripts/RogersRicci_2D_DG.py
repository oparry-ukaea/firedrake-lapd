# RogersRicci_2D_DG.py

from firedrake import *
import math
from irksome import Dt, MeshConstant, TimeStepper, GaussLegendre
import time
import numpy

meshres = 128
#mesh = SquareMesh(meshres, meshres, 100, quadrilateral=True)
mesh = Mesh("test_circle.msh")
V1 = FunctionSpace(mesh, "DG", 1) # for vorticity (w)
V2 = FunctionSpace(mesh, "DG", 1) # for temperature (T)
V3 = FunctionSpace(mesh, "CG", 1) # for potential (phi)
V4 = VectorFunctionSpace(mesh, "CG", 1)  # for drift velocity
V = V1*V2

T_end = 0.50
timeres = 500

t = Constant(0.0)
dt = Constant(T_end/timeres)

# parameters for irksome
butcher_tableau = GaussLegendre(1)

# model parameters (TODO put them here)

x, y = SpatialCoordinate(mesh)
wT = Function(V)
w, T = split(wT)
v1, v2 = TestFunctions(V)
phi = TrialFunction(V3)
v3 = TestFunction(V3)

driftvel_out = Function(V4)

Lphi = inner(grad(phi),grad(v3))*dx
Rphi = -w*v3*dx
phi_s = Function(V3)
bc1 = DirichletBC(V3, 0, 'on_boundary')

ST = Function(wT.sub(1))  # source T
rs = 1.0*0.28*100/1.4  # 0.28 in paper, lengthscale l_0 is 0.014m = rho_{s0} in the paper
Ls = 5.0*0.007*100/1.4 #TRIALCODE uplift (first factor) can be used to decrease the gradients
ST.interpolate(0.5*(1-tanh((sqrt((x-50.0)*(x-50.0)+(y-50.0)*(y-50.0))-rs)/Ls)))

# TRIALCODE interpolate the steady-state solution
# phi is determined from this so does not need to be initialized
eps_init = 1.0e-6

def r(x,y):
  return(sqrt((x-50)**2+(y-50)**2+eps_init))

def farg(x,y):
   return((sqrt((x-50)**2+(y-50)**2)-rs)/Ls)

wT.sub(0).interpolate(1.5*(38.2/0.0278)*((2*tanh(farg(x,y))/(Ls**2*cosh(farg(x,y))**2))-(1/(Ls*r(x,y)*cosh(farg(x,y))**2))))
wT.sub(1).interpolate((38.2/0.0278)*0.5*(1-tanh(farg(x,y))))

testval = assemble((wT.sub(0))*dx)  # ought to be zero
print("integral of vorticity over domain is "+str(float(testval))+"\n")

driftvel = as_vector([2.5*grad(phi_s)[1],-2.5*grad(phi_s)[0]])

# TRIALCODE check init data
#File("RogersRicci_2D_DG.pvd").write(wT.sub(0), wT.sub(1))
#quit()


# eps is needed to stop immediate crash
eps = 0.01

# needed to stop div zero
T_eps = 0.01  # was 0.01

norm = FacetNormal(mesh)
driftvel_n = 0.5*(dot(driftvel, norm)+abs(dot(driftvel, norm)))

#TRIALCODE exp(y) replaced with 1+y
# note should probably remove the 0.00*exp(stuff), hopefully the code generator was smart enough to remove them when I ran this ...
F = Dt(w)*v1*dx + Dt(T)*v2*dx \
  - v1*div(w*driftvel)*dx - v2*div(T*driftvel)*dx\
  - 0.00667*(1-(1+3-phi_s/sqrt(T*T+T_eps*T_eps))/(1+0.00*exp(3-phi_s/sqrt(T*T+T_eps*T_eps))))*(v1)*dx \
  + (0.0278*T*(1.71*(1+3-phi_s/sqrt(T*T+T_eps*T_eps))/(1+0.00*exp(3-phi_s/sqrt(T*T+T_eps*T_eps)))-0.71)-38.2*ST)*(v2)*dx \
  + driftvel_n('-')*(w('-') - w('+'))*v1('-')*dS \
  + driftvel_n('+')*(w('+') - w('-'))*v1('+')*dS \
  + driftvel_n('-')*(T('-') - T('+'))*v2('-')*dS \
  + driftvel_n('+')*(T('+') - T('-'))*v2('+')*dS \

# params taken from Cahn-Hilliard example cited above
params = {'snes_monitor': None, 'snes_max_it': 100,
          'snes_linesearch_type': 'l2',
          'ksp_type': 'preonly',
          'pc_type': 'lu', 'mat_type': 'aij',
          'pc_factor_mat_solver_type': 'mumps'}

stepper = TimeStepper(F, butcher_tableau, t, dt, wT, solver_parameters=params)

# this is intended to be direct solver
linparams = {"mat_type": "aij",
          "snes_type": "ksponly",
          "ksp_type": "preonly",
          "pc_type": "lu"}

outfile = File("RogersRicci_2D_DG.pvd")

cnt=0
start = time.time()

while float(t) < float(T_end):
    if (float(t) + float(dt)) >= T_end:
        dt.assign(T_end - float(t))
    solve(Lphi==Rphi, phi_s, solver_parameters=linparams, bcs=bc1)
    print("solved phi\n")
    driftvel = as_vector([2.5*grad(phi_s)[1],-2.5*grad(phi_s)[0]])
    driftvel_n = 0.5*(dot(driftvel, norm)+abs(dot(driftvel, norm)))
    driftvel_out.interpolate(driftvel)  # TRIALCODE for output driftvel
    stepper.advance()
    t.assign(float(t) + float(dt))
    print(float(t), float(dt))
    testval = assemble((wT.sub(0))*dx)
    print("integral of vorticity over domain is "+str(float(testval))+"\n")

    if(cnt % 10 == 0):
       print("outputting data ...\n")
       ws,Ts = wT.split()
       ws.rename("vorticity")
       Ts.rename("temperature")
       phi_s.rename("potential")
       driftvel_out.rename("driftvel")
       outfile.write(ws, Ts, phi_s, driftvel_out)
    cnt=cnt+1


end = time.time()
wall_time = end-start

print("done.")
print("\n")
print("wall time:"+str(wall_time)+"\n")
