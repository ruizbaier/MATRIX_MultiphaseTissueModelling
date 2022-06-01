'''

Multiphase (mixture-theory based) description of a tumour model

dc is the cells' displacement
alpha is the Lagrangian solid volume fraction
pw is the interstitial fluid pressure
eta is the oxygen concentration 

the geometry is a glioma at a given stage

Boundary conditions are of pure-stress type. Therefore an additional step is needed for sake of uniqueness of solutions 

the equations in strong form read as follows

(balance of momentum for the mixture) 
- Div (P) = 0, with P = partial(W)/partial(F) - pw*J*F^{-T} - alpha*Sigma(alpha)*J*F^{-T}

(balance of mass for the solid)
d/dt(J*alpha) = J * qc(alpha,eta)

(balance of mass for the mixture)
d/dt(J) - Div(J/k*F^{-1}*F^{-T}*Grad(pw)) = 0 

(oxygen diffusion)
- Div(J*F^{-1}*F^{-T}*Grad(eta)) = J*f(alpha,eta)

'''


from fenics import *

parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["allow_extrapolation"]= True
parameters["form_compiler"]["cpp_optimize"] = True

Sigma  = lambda alpha: Heavi(alpha-alphamin)*(alpha-alphastar)*pow(1-alpha,-2)
qc     = lambda alpha,C: alpha*(1-alpha)*(1+s1)*C/(1+s1*C) - (s2+s3*C)*alpha/(1+s4*C)
f      = lambda alpha,C: -Q*alpha*C/(1+Qhat*C)
Jc_fun = lambda d: det(I + grad(d))

Heavi  = lambda u: conditional(ge(u, 0.0), 1, 0)

# ******* Constants
ndim     = 2
alphamin = Constant(0.8)
alphastar= Constant(0.8)
alpha0   = Constant(0.8)
mu       = Constant(1.)
lmbda    = Constant(1.)
pwout    = Constant(0.)
Cout     = Constant(1.)
k        = Constant(1.)
zero     = Constant((0,0))
Q        = Constant(0.5)
Qhat     = Constant(0.)
s1       = Constant(10.)
s2       = Constant(0.5)
s3       = Constant(0.5)
s4       = Constant(10.)

# ******* Define mesh and define function spaces ****** #
mesh = Mesh("meshes/glioma.xml")
N = FacetNormal(mesh)

dx = Measure("dx", domain=mesh)
#ds = Measure("ds", subdomain_data = bdry)

#output files and options 
fileO = XDMFFile(mesh.mpi_comm(), "outputs/MultiphaseGlioma.xdmf")
fileO.parameters['rewrite_function_mesh']=False
fileO.parameters["functions_share_mesh"] = True
fileO.parameters["flush_output"] = True

# ********* Time constants ********* #

t = 0.0; dt = 0.5; Tf = 100.; freqSave = 20

# ********* Finite dimensional spaces ********* #

P1 = FiniteElement("CG", mesh.ufl_cell(), 1)
Bub  = FiniteElement("Bubble", mesh.ufl_cell(), 3)
P1b  = VectorElement(P1 + Bub)
RM   = VectorElement('R', triangle, 0, dim=3)
#P2v = VectorElement("CG", mesh.ufl_cell(), 2)
Hh = FunctionSpace(mesh,MixedElement([P1b, P1, P1, P1, RM]))


nullspace=[Constant((1,0)), Constant((0,1)),\
               Expression(('-x[1]','x[0]'),degree = 1)]


print("**************** Total DoFs = ", Hh.dim())

# trial and test functions

Sol  = Function(Hh); dSol = TrialFunction(Hh); dTest = TestFunction(Hh)
dc,           alpha,      pw,      C, s_chi = split(Sol)
dc_test, alpha_test, pw_test, C_test, s_xi  = split(dTest)

# ********* Initial conditions ******* #

dc_old = interpolate(project(zero,Hh.sub(0).collapse()),Hh.sub(0).collapse())
alpha_old = interpolate(alpha0,Hh.sub(1).collapse())

# ********* Boundaries and boundary conditions ******* #
# Pure traction for the mechanics, pure Dirichlet for pressure and oxygen
bcP  = DirichletBC(Hh.sub(2), pwout, 'on_boundary')
bcE  = DirichletBC(Hh.sub(3), Cout, 'on_boundary')
bcs = [bcP,bcE]

# ******** Define hyperelastic quantities  ************* #

I = Identity(ndim); Fc = I + grad(dc); Fc = variable(Fc)
Bc = Fc*Fc.T; Cc = Fc.T*Fc; Jc = det(Fc); invFc = inv(Fc); I1 = tr(Cc)

# Neo-Hookean strain energy density
Wc = 0.5*mu*(I1-ndim-2*ln(Jc)) + 0.5*lmbda*(Jc-1)**2

# Effective 1st Piola-Kirchhoff stress tensor for the cell phase
Pc = diff(Wc,Fc)

# total stress (assuming an active stress approach)
P = Pc - (pw + alpha*Sigma(alpha)) * Jc * invFc.T

# ********  Weak form ********** #
     
FF = inner(P, grad(dc_test)) * dx \
     + ((Jc-Jc_fun(dc_old))*alpha+(alpha-alpha_old)*Jc)/dt * alpha_test * dx \
     - Jc*qc(alpha,C) * alpha_test * dx \
     + (Jc-Jc_fun(dc_old))/dt*pw_test * dx \
     + dot(Jc/k*inv(Cc)*grad(pw),grad(pw_test)) * dx \
     + dot(Jc*inv(Cc)*grad(C),grad(C_test))*dx \
     - Jc*f(alpha,C)*C_test * dx

for i, ns_i in enumerate(nullspace):
    chi = s_chi[i]
    xi  = s_xi[i]
    FF += chi*inner(dc_test, ns_i)*dx + xi*inner(dc, ns_i)*dx



Tang = derivative(FF, Sol, dSol)
problem = NonlinearVariationalProblem(FF, Sol, bcs, J=Tang)
solver  = NonlinearVariationalSolver(problem)
solver.parameters['nonlinear_solver']                    = 'newton'
solver.parameters['newton_solver']['linear_solver']      = 'mumps'
solver.parameters['newton_solver']['absolute_tolerance'] = 1e-6
solver.parameters['newton_solver']['relative_tolerance'] = 1e-6
solver.parameters['newton_solver']['maximum_iterations'] = 12

# ************* Time loop ********** #
inc = 0;

while (t <= Tf):

    print("t=%.3f" % t)

    solver.solve()
    dch,alphah,pwh,Ch, chih = Sol.split()
    assign(alpha_old,alphah)
    assign(dc_old,dch)
    
    if (inc % freqSave == 0):
        
        dch.rename("dc","dc"); fileO.write(dch,t)
        alphah.rename("alpha","alpha"); fileO.write(alphah,t)
        pwh.rename("pw","pw"); fileO.write(pwh,t)
        Ch.rename("C","C"); fileO.write(Ch,t)

    t += dt; inc += 1

# ************* End **************** #    
