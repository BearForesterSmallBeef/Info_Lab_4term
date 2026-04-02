from basix.ufl import element
from ufl import (
    Coefficient, Constant, FunctionSpace, Mesh, TestFunction, TrialFunction,
    Identity, dx, ds, grad, inner, sym, div, dot, FacetNormal
)

# 1. Сетка
coord_element = element("Lagrange", "tetrahedron", 1, shape=(3,))
mesh = Mesh(coord_element)

# 2. Пространства
e_v = element("Lagrange", "tetrahedron", 2, shape=(3,))
e_q = element("Lagrange", "tetrahedron", 1)

V = FunctionSpace(mesh, e_v)
Q = FunctionSpace(mesh, e_q)

# Искомые функции для левых частей (Trial)
u = TrialFunction(V)
p = TrialFunction(Q)

# Тестовые функции
v = TestFunction(V)
q = TestFunction(Q)

# Коэффициенты для правых частей (Known)
u_n = Coefficient(V)      # Скорость на шаге n
u_s = Coefficient(V)      # Предварительная скорость u* (результат Шага 1)
p_n = Coefficient(Q)      # Давление на шаге n
p_s = Coefficient(Q)      # Давление на шаге n+1 (результат Шага 2)

dt = Constant(mesh)
rho = Constant(mesh)      # Плотность
mu = Constant(mesh)       # Динамическая вязкость

# =================================================================
# ШАГ 1: Предварительная скорость u*
# =================================================================
a1 = (rho / dt) * inner(u, v) * dx \
   + rho * inner(dot(grad(u), u_n), v) * dx \
   + mu * inner(grad(u), grad(v)) * dx

L1 = (rho / dt) * inner(u_n, v) * dx \
   - inner(grad(p_n), v) * dx

# =================================================================
# ШАГ 2: Давление p_next (Уравнение Пуассона)
# Добавлен grad(p_n) для инкрементальной схемы (IPCS)
# =================================================================
a2 = inner(grad(p), grad(q)) * dx
L2 = inner(grad(p_n), grad(q)) * dx - (rho / dt) * inner(div(u_s), q) * dx
# =================================================================
# ШАГ 3: Коррекция скорости u_next
# =================================================================
a3 = inner(u, v) * dx
# ИСПРАВЛЕНИЕ: Используем p_s (коэффициент), а не p (TrialFunction)
L3 = inner(u_s, v) * dx - (dt / rho) * inner(grad(p_s - p_n), v) * dx

# =================================================================
# ЛОБОВОЕ СОПРОТИВЛЕНИЕ (Drag Force)
# =================================================================
def epsilon(u_vec): return sym(grad(u_vec))
def sigma(u_vec, p_scal): return 2.0 * mu * epsilon(u_vec) - p_scal * Identity(len(u_vec))

n = FacetNormal(mesh)
traction = dot(sigma(u_n, p_n), n)
drag_force = traction[0] * ds(4)

# Итоговый список форм для C++
forms = [a1, L1, a2, L2, a3, L3, drag_force]