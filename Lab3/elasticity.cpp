#include "elasticity.h"
#include <basix/finite-element.h>
#include <cmath>
#include <dolfinx.h>
#include <dolfinx/fem/Constant.h>
#include <dolfinx/fem/CoordinateElement.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/io/XDMFFile.h>
#include <dolfinx/la/petsc.h>
#include <petscmat.h>
#include <petscsys.h>
#include <petscsystypes.h>
#include <utility>
#include <vector>
#include <map>
#include <iostream>
#include <span>
#include <optional>

using namespace dolfinx;
using T = PetscScalar;
using U = typename dolfinx::scalar_value_t<T>;

int main(int argc, char* argv[])
{
  dolfinx::init_logging(argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);

  {
    // load om_nom network
    fem::CoordinateElement<U> coord_element(mesh::CellType::tetrahedron, 1);
    
    io::XDMFFile file(MPI_COMM_WORLD, "om_nom.xdmf", "r");
    auto mesh = std::make_shared<mesh::Mesh<U>>(
        file.read_mesh(coord_element, mesh::GhostMode::shared_facet, "Grid")
    );
    file.close();

    // векторное пространство
    auto basix_element = basix::create_element<U>(
        basix::element::family::P, basix::cell::type::tetrahedron, 1,
        basix::element::lagrange_variant::unset,
        basix::element::dpc_variant::unset, false);

    // Указываем векторную структуру {3}
    std::vector<std::size_t> value_shape = {3};
    auto element = std::make_shared<fem::FiniteElement<U>>(basix_element, value_shape, false);

    auto V = std::make_shared<fem::FunctionSpace<U>>(fem::create_functionspace<U>(mesh, element));

    // константы моделирования
    auto mu = std::make_shared<fem::Constant<T>>(5000.0);
    auto lmbda = std::make_shared<fem::Constant<T>>(5000.0);
    auto f = std::make_shared<fem::Function<T>>(V);

    // суета с типами
    auto V_const = std::shared_ptr<const fem::FunctionSpace<U>>(V);
    auto mu_const = std::shared_ptr<const fem::Constant<T>>(mu);
    auto lmbda_const = std::shared_ptr<const fem::Constant<T>>(lmbda);
    auto f_const = std::shared_ptr<const fem::Function<T>>(f);

    std::vector<std::shared_ptr<const fem::FunctionSpace<U>>> spaces_a = {V_const, V_const};
    std::map<std::string, std::shared_ptr<const fem::Constant<T>>> consts_a = {
        {"mu", mu_const}, 
        {"lmbda", lmbda_const}
    };
    
    std::vector<std::shared_ptr<const fem::FunctionSpace<U>>> spaces_L = {V_const};
    std::map<std::string, std::shared_ptr<const fem::Function<T>>> coeffs_L = {
        {"f", f_const}
    };

    fem::Form<T> a = fem::create_form<T>(*form_elasticity_a, spaces_a, {}, consts_a, {}, {});
    fem::Form<T> L = fem::create_form<T>(*form_elasticity_L, spaces_L, coeffs_L, {}, {}, {});

    // граничные условия
    auto x_span = mesh->geometry().x(); 
    std::size_t num_vertices = x_span.size() / 3; 
    
    U min_z = 1e9;
    for (std::size_t i = 0; i < num_vertices; ++i) {
        if (x_span[3*i + 2] < min_z) min_z = x_span[3*i + 2];
    }
    std::cout << "Min Z found: " << min_z << std::endl;

    auto ground_check = [min_z](auto x) {
          std::vector<std::int8_t> marker(x.extent(1), false);
          for (std::size_t p = 0; p < x.extent(1); ++p) {
            if (std::abs(x(2, p) - min_z) < 0.05) marker[p] = true; 
          }
          return marker;
    };

    std::vector facets = mesh::locate_entities_boundary(*mesh, 2, ground_check);
    std::vector bdofs = fem::locate_dofs_topological(*V->mesh()->topology_mutable(), *V->dofmap(), 2, facets);
    
    // чото от нейронки без этого не робит
    auto bc_func = std::make_shared<fem::Function<T>>(V);
    bc_func->interpolate([](auto x) -> std::pair<std::vector<T>, std::vector<std::size_t>> {
        std::size_t N = x.extent(1);
        return {std::vector<T>(3 * N, 0.0), {3, N}};
    });

    // мотеша тун тун
    fem::DirichletBC<T> bc(std::shared_ptr<const fem::Function<T>>(bc_func), bdofs);

    la::petsc::Matrix A(fem::petsc::create_matrix(a), false);

    MatZeroEntries(A.mat());
    fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES), a, {bc});
    MatAssemblyBegin(A.mat(), MAT_FLUSH_ASSEMBLY);
    MatAssemblyEnd(A.mat(), MAT_FLUSH_ASSEMBLY);
    fem::set_diagonal<T>(la::petsc::Matrix::set_fn(A.mat(), INSERT_VALUES), *V, {bc});
    MatAssemblyBegin(A.mat(), MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A.mat(), MAT_FINAL_ASSEMBLY);

    auto u = std::make_shared<fem::Function<T>>(V);
    la::Vector<T> b(L.function_spaces()[0]->dofmap()->index_map,
                    L.function_spaces()[0]->dofmap()->index_map_bs());

    la::petsc::KrylovSolver lu(MPI_COMM_WORLD);
    la::petsc::options::set("ksp_type", "preonly");
    la::petsc::options::set("pc_type", "lu");
    lu.set_from_options();

    lu.set_operator(A.mat());

    io::VTKFile file_out(MPI_COMM_WORLD, "om_nom_animation.pvd", "w");

    int steps = 240;              // Количество кадров
    double max_gravity = -200.0;

    for (int step = 0; step <= steps; ++step) {
        double t = static_cast<double>(step) / steps;
        double current_gravity = max_gravity * t;

        std::cout << "Step: " << step << " / " << steps << ", Gravity: " << current_gravity << std::endl;

        // 1. Обновляем нагрузку (f)
        f->interpolate([current_gravity](auto x) -> std::pair<std::vector<T>, std::vector<std::size_t>> {
              std::size_t N = x.extent(1);
              std::vector<T> f_vals(3 * N);
              for (std::size_t p = 0; p < N; ++p) {
                f_vals[0 * N + p] = 0.0;
                f_vals[1 * N + p] = 0.0;
                f_vals[2 * N + p] = current_gravity;
              }
              return {f_vals, {3, N}};
        });

        // 2. Пересобираем вектор b (так как f изменилась)
        std::ranges::fill(b.array(), 0);
        fem::assemble_vector(b.array(), L);
        
        // 3. Применяем граничные условия к вектору b
        fem::apply_lifting(b.array(), {a}, {{bc}}, {}, T(1));
        b.scatter_rev(std::plus<T>());
        bc.set(b.array(), std::nullopt);

        // 4. Решаем
        la::petsc::Vector _u(la::petsc::create_vector_wrap(*u->x()), false);
        la::petsc::Vector _b(la::petsc::create_vector_wrap(b), false);
        lu.solve(_u.vec(), _b.vec());

        // 5. Обновляем ghost-значения и сохраняем кадр
        u->x()->scatter_fwd();
        
        file_out.write(std::vector<std::reference_wrapper<const fem::Function<T>>>{*u}, t);
    }
  }
   

  PetscFinalize();
  return 0;
}