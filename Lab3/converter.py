import meshio
import numpy as np

# Читаем твою сетку
msh = meshio.read("stl.msh")

# Нам нужны только тетраэдры (объем) для физики
tets = msh.get_cells_type("tetra")

# Создаем новую сетку только из тетраэдров
mesh_out = meshio.Mesh(points=msh.points, cells={"tetra": tets})

# Сохраняем в XDMF (формат, который любит FEniCS)
meshio.write("om_nom.xdmf", mesh_out)

print("Сетка успешно сконвертирована в om_nom.xdmf!")