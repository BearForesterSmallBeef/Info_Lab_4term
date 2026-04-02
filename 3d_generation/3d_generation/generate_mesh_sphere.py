import gmsh
import sys


def generate_channel_with_sphere():
    gmsh.initialize()
    gmsh.model.add("ChannelWithSphere")

    occ = gmsh.model.occ

    L, W, H = 2.0, 0.6, 0.6  # Длина, Ширина, Высота
    # Параметры сферы
    R = 0.1
    cx, cy, cz = 0.7, 0.3, 0.3 # Координаты центра сферы

    channel = occ.addBox(0, 0, 0, L, W, H)
    sphere = occ.addSphere(cx, cy, cz, R)

    # вычиатем сферу из трубы
    # первый аргумент - объекты "из которых" вычитаем
    # второй - объекты "которые" вычитаем
    # dimTags - это список кортежей [(размерность, id)]
    fluid_domain, _ = occ.cut([(3, channel)], [(3, sphere)])

    occ.synchronize()

    # Нам нужно найти поверхности, чтобы потом задать на них условия
    # Мы ищем их по координатам их центров

    # ищем вход - грань при x = 0
    inlet_list = gmsh.model.getEntitiesInBoundingBox(-0.01, -0.01, -0.01, 0.01, W + 0.01, H + 0.01, 2)
    # ищем выход - грань при x = L
    outlet_list = gmsh.model.getEntitiesInBoundingBox(L - 0.01, -0.01, -0.01, L + 0.01, W + 0.01, H + 0.01, 2)
    #ищем сферу ящиком
    sphere_list = gmsh.model.getEntitiesInBoundingBox(cx - R - 0.01, cy - R - 0.01, cz - R - 0.01,
                                                      cx + R + 0.01, cy + R + 0.01, cz + R + 0.01, 2)
    inlet_ids = [tag for dim, tag in inlet_list]
    outlet_ids = [tag for dim, tag in outlet_list]
    sphere_surface = [tag for dim, tag in sphere_list]

    # поверхность сферы - это то, что осталось после вычитания, но не вход, не выход и не внешние стенки
    # найдем сферу как поверхность, чьи координаты центра близки к (cx, cy, cz)
    all_surfaces = gmsh.model.getEntities(2)
    wall_surfaces = []

    for dim, tag in all_surfaces:
        if tag not in inlet_ids and tag not in outlet_ids and tag not in sphere_surface:
            wall_surfaces.append(tag)

    gmsh.model.addPhysicalGroup(2, inlet_ids, 1, name="Inlet")
    gmsh.model.addPhysicalGroup(2, outlet_ids, 2, name="Outlet")
    gmsh.model.addPhysicalGroup(2, wall_surfaces, 3, name="Walls")
    gmsh.model.addPhysicalGroup(2, sphere_surface, 4, name="Sphere")

    # объем самого жижи
    gmsh.model.addPhysicalGroup(3, [fluid_domain[0][1]], 10, name="Fluid")

    print(f"Inlet IDs: {inlet_ids}")
    print(f"Outlet IDs: {outlet_ids}")
    print(f"Sphere IDs: {sphere_surface}")
    print(f"Wall IDs: {wall_surfaces}")

    base_size_min = 0.01
    base_size_max = 0.08
    # base_size_min = 0.005
    # base_size_max = 0.04
    gmsh.option.setNumber("Mesh.CharacteristicLengthMin", base_size_min)
    gmsh.option.setNumber("Mesh.CharacteristicLengthMax", base_size_max)

    dist_field = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist_field, "FacesList", sphere_surface)

    threshold_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold_field, "InField", dist_field)

    gmsh.model.mesh.field.setNumber(threshold_field, "InField", dist_field)
    # Размер ячейки у поверхности сферы
    gmsh.model.mesh.field.setNumber(threshold_field, "SizeMin", base_size_min)
    # Размер ячейки вдали от сферы
    gmsh.model.mesh.field.setNumber(threshold_field, "SizeMax", base_size_max)
    # Расстояние, до которого действует SizeMin
    # gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", 0.025)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", 0.05)
    # Расстояние, на котором размер достигает SizeMax
    # gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", 0.20)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", 0.20)

    gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

    gmsh.model.mesh.generate(3)
    gmsh.write("mesh/channel_sphere.msh")
    gmsh.finalize()


if __name__ == "__main__":
    generate_channel_with_sphere()
