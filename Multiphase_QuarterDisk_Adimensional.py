'''


multiphase (mixture-theory based) description of a tumour model

dc is the cells' displacement
alpha is the Lagrangian solid volume fraction
pw is the interstitial fluid pressure
eta is the oxygen concentration 

the geometry is simply a quarter of a disk

the equations in strong form read as follows


- Div (P) = 0, with P = partial(W)/partial(F) - pw*J*F^{-T} - alpha*Sigma(alpha)*J*F^{-T}
d/dt(J) - Div(J*F^{-1}*F^{-T}*(1-alpha)/(k1*alpha)*Grad(pw)) = 0 
d/dt(J*alpha) = J * qc(alpha,eta) 
- Div(J*F^{-1}*F^{-T}*Grad(eta)) = J*f(alpha,eta)


boundary conditions are of sliding type on the straight segments and of zero-traction on the circular arc
'''


from fenics import *

parameters["form_compiler"]["representation"] = "uflacs"
parameters["form_compiler"]["quadrature_degree"] = 4
parameters["allow_extrapolation"]= True
parameters["form_compiler"]["cpp_optimize"] = True

Sigma  = lambda alpha: Heavi(alpha-alphamin)*gamma*(alpha-alphastar)*pow(1-alpha,-2)
qc     = lambda alpha,C: alpha*(1-alpha)*(1+s1)*C/(1+s1*C) - (s2+s3*C)*alpha/(1+s4*C)
f      = lambda alpha,C: -Q*alpha*C/(1+hatQ1*C)
Jc_fun = lambda d: det(I + grad(d))

Heavi  = lambda u: conditional(ge(u, 0.0), 1, 0)

# ******* Constants
ndim     = 2
alphamin = Constant(0.8)
alphastar= Constant(0.8)
alpha0   = Constant(0.5)
pwout    = Constant(0.)
Cout     = Constant(1.)
k        = Constant(0.25)
zero     = Constant((0,0))
Q        = Constant(0.5)
hatQ1    = Constant(0.)
gamma    = Constant(1.0)
s1       = Constant(10.)
s2       = Constant(0.5)
s3       = Constant(0.5)
s4       = Constant(10.)

mu       = Constant(0.25)
lmbda    = Constant(1.)

# Young modulus, Poisson ratio
#E, nu = 1.e2, 0.45 # N/m^2, [-]
#mu,lmbda = Constant(E/(2*(1+nu))), Constant(E*nu/((1+nu)*(1-2*nu))) # N/m^2

# ******* Define mesh and define function spaces ****** #
mesh = Mesh("meshes/quarterDisk.xml")
bdry = bdry = MeshFunction("size_t", mesh, "meshes/quarterDisk_facet_region.xml")
circ =32; left= 33; bot = 31
N = FacetNormal(mesh)

dx = Measure("dx", domain=mesh)
ds = Measure("ds", subdomain_data = bdry)

#output files and options 
fileO = XDMFFile(mesh.mpi_comm(), "outputs/MultiphaseQuarterDisk.xdmf")
fileO.parameters['rewrite_function_mesh']=False
fileO.parameters["functions_share_mesh"] = True
fileO.parameters["flush_output"] = True

# ********* Time constants ********* #

t = 0.0; dt = 0.5; Tf = 40.; freqSave = 1

# ********* Finite dimensional spaces ********* #

P1 = FiniteElement("CG", mesh.ufl_cell(), 1)
Bub  = FiniteElement("Bubble", mesh.ufl_cell(), 3)
P1b  = VectorElement(P1 + Bub)
#P2v = VectorElement("CG", mesh.ufl_cell(), 2)
Hh = FunctionSpace(mesh,MixedElement([P1b, P1, P1, P1]))

print("**************** Total DoFs = ", Hh.dim())

# trial and test functions

Sol  = Function(Hh); dSol = TrialFunction(Hh); dTest = TestFunction(Hh)
dc,           alpha,      pw,      C = split(Sol)
dc_test, alpha_test, pw_test, C_test = split(dTest)

# ********* Initial conditions ******* #

dc_old = interpolate(project(zero,Hh.sub(0).collapse()),Hh.sub(0).collapse())
alpha_old = interpolate(alpha0,Hh.sub(1).collapse())

# ********* Boundaries and boundary conditions ******* #

bcU1 = DirichletBC(Hh.sub(0).sub(1), project(Constant(0.),Hh.sub(0).sub(1).collapse()), bdry, bot)
bcU2 = DirichletBC(Hh.sub(0).sub(0), project(Constant(0.),Hh.sub(0).sub(0).collapse()), bdry, left)
bcP  = DirichletBC(Hh.sub(2), pwout, bdry, circ)
bcE  = DirichletBC(Hh.sub(3), Cout, bdry, circ)

bcs = [bcU1,bcU2,bcP,bcE]

isBlocked = Expression("x[1] < x[0]? 1:0", degree = 0) 

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
     + dot(Jc*(k*alpha)/(1-alpha)*inv(Cc)*grad(pw),grad(pw_test)) * dx \
     + dot(Jc*inv(Cc)*grad(C),grad(C_test))*dx \
     - Jc*f(alpha,C)*C_test * dx 

Tang = derivative(FF, Sol, dSol)
problem = NonlinearVariationalProblem(FF, Sol, bcs, J=Tang)
solver  = NonlinearVariationalSolver(problem)
solver.parameters['nonlinear_solver']                    = 'newton'#snes
solver.parameters['newton_solver']['linear_solver']      = 'petsc'#mumps
solver.parameters['newton_solver']['absolute_tolerance'] = 1e-6
solver.parameters['newton_solver']['relative_tolerance'] = 1e-6
solver.parameters['newton_solver']['maximum_iterations'] = 15

# ************* Time loop ********** #
inc = 0;

while (t <= Tf):

    print("t=%.3f" % t)

    solver.solve()
    dch,alphah,pwh,Ch = Sol.split()
    assign(alpha_old,alphah)
    assign(dc_old,dch)
    
    if (inc % freqSave == 0):
        
        dch.rename("dc","dc"); fileO.write(dch,t)
        alphah.rename("alpha","alpha"); fileO.write(alphah,t)
        pwh.rename("pw","pw"); fileO.write(pwh,t)
        Ch.rename("C","C"); fileO.write(Ch,t)

    t += dt; inc += 1

# ************* End **************** #    
