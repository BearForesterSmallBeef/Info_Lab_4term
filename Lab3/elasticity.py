from basix.ufl import element
from ufl import (
    Coefficient,
    Constant,
    FunctionSpace,
    Mesh,
    TestFunction,
    TrialFunction,
    Identity,
    dx,
    grad,
    inner,
    sym,
    tr,
)

coord_element = element("Lagrange", "tetrahedron", 1, shape=(3,))
mesh = Mesh(coord_element)

e = element("Lagrange", "tetrahedron", 1, shape=(3,))
V = FunctionSpace(mesh, e)

u = TrialFunction(V)
v = TestFunction(V)
f = Coefficient(V)

mu = Constant(mesh)
lmbda = Constant(mesh)

def epsilon(u):
    return sym(grad(u))

def sigma(u):
    return 2.0 * mu * epsilon(u) + lmbda * tr(epsilon(u)) * Identity(len(u))

a = inner(sigma(u), epsilon(v)) * dx
L = inner(f, v) * dx
