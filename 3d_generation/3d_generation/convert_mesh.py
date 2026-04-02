import meshio


"""
def convert():
    msh = meshio.read("mesh/channel_sphere.msh")

    # Сохраняем только тетраэдры для объема
    tetra_mesh = meshio.Mesh(points=msh.points,
                             cells={"tetra": msh.cells_dict["tetra"]},
                             cell_data={"f": [msh.cell_data_dict["gmsh:physical"]["tetra"]]})
    meshio.write("mesh/sphere_volume.xdmf", tetra_mesh)

    # Сохраняем треугольники для границ (чтобы C++ видел Inlet, Outlet и т.д.)
    triangle_mesh = meshio.Mesh(points=msh.points,
                                cells={"triangle": msh.cells_dict["triangle"]},
                                cell_data={"f": [msh.cell_data_dict["gmsh:physical"]["triangle"]]})
    meshio.write("mesh/sphere_facets.xdmf", triangle_mesh)
"""


def convert():
    msh = meshio.read("mesh/channel_sphere.msh")

    # Для объема
    tetra_mesh = meshio.Mesh(points=msh.points,
                             cells={"tetra": msh.cells_dict["tetra"]},
                             cell_data={"Grid": [msh.cell_data_dict["gmsh:physical"]["tetra"]]}) # Тут Grid
    meshio.write("mesh/sphere_volume.xdmf", tetra_mesh)

    # Для границ
    triangle_mesh = meshio.Mesh(points=msh.points,
                                cells={"triangle": msh.cells_dict["triangle"]},
                                cell_data={"Grid": [msh.cell_data_dict["gmsh:physical"]["triangle"]]}) # Тут Grid
    meshio.write("mesh/sphere_facets.xdmf", triangle_mesh)


if __name__ == "__main__":
    convert()
