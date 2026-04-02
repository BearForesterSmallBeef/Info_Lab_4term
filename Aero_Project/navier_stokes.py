from basix.ufl import element, mixed_element
from ufl import (
    Coefficient, Constant, FunctionSpace, Mesh, TestFunctions,
    split, Identity, dx, grad, inner, sym, div, dot, derivative, ds, FacetNormal
)

# сетка (Тетраэдры, shape=(3,))
coord_element = element("Lagrange", "tetrahedron", 1, shape=(3,))
mesh = Mesh(coord_element)

e_u = element("Lagrange", "tetrahedron", 2, shape=(3,))
e_p = element("Lagrange", "tetrahedron", 1)
e_W = mixed_element([e_u, e_p])
W = FunctionSpace(mesh, e_W)

w = Coefficient(W)
u, p = split(w)
v, q = TestFunctions(W)

nu = Constant(mesh)


def epsilon(u):
    return sym(grad(u))


def sigma(u, p):
    return 2.0 * nu * epsilon(u) - p * Identity(len(u))


# вектор невязки F
F = (inner(dot(u, grad(u)), v) * dx
   + inner(sigma(u, p), epsilon(v)) * dx
   + q * div(u) * dx)

# якобиан (Производная F по w)
J = derivative(F, w)

# интергал силы сопротивления (Drag Force)
n = FacetNormal(mesh)
# Вектор поверхностных сил (Traction vector) T = sigma * n
traction = dot(sigma(u, p), n)

# Лобовое сопротивление - это проекция вектора сил на ось X (индекс 0).
# Интегрируем по поверхности шарика (маркер 3)
drag_force = traction[0] * ds(3)
# Указываем ffcx, что хотим скомпилировать эти формы в плюсы
forms =[F, J, drag_force]