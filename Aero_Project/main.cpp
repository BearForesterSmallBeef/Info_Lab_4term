#include "navier_stokes_chorin.h"
#include <basix/finite-element.h>
#include <dolfinx.h>
#include <dolfinx/fem/Constant.h>
#include <dolfinx/fem/CoordinateElement.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/la/petsc.h>
#include <dolfinx/io/XDMFFile.h>
#include <petscmat.h>
#include <petscsys.h>
#include <iostream>
#include <span>
#include <vector>
#include <map>
#include <algorithm>
#include <optional>
#include <fstream>

using namespace dolfinx;
using T = PetscScalar;
using U = typename dolfinx::scalar_value_t<T>;

int main(int argc, char* argv[])
{
    dolfinx::init_logging(argc, argv);
    PetscInitialize(&argc, &argv, nullptr, nullptr);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::ofstream log_file;
    std::streambuf* original_cout_buf = std::cout.rdbuf(); // сохраняем оригинальный вывод в консоль

    if (rank == 0) {
        log_file.open("log.txt", std::ios::out);
        // перенаправляем std::cout в наш файл log.txt
        std::cout.rdbuf(log_file.rdbuf()); 
        std::cout << "=== НАЧАЛО РАСЧЕТА CFD ===" << std::endl;
    }

    {
        fem::CoordinateElement<U> cmap(mesh::CellType::tetrahedron, 1);

        // читаем сетку
        io::XDMFFile file_vol(MPI_COMM_WORLD, "mesh/sphere_volume.xdmf", "r");
        // XPath ("/Xdmf/Domain")
        auto mesh = std::make_shared<mesh::Mesh<U>>(
            file_vol.read_mesh(cmap, mesh::GhostMode::shared_facet, "Grid")
        );
        file_vol.close();
        mesh->topology()->create_entities(2); 
        mesh->topology()->create_connectivity(2, 3);

        // читаем теги сферы
        io::XDMFFile file_facets(MPI_COMM_WORLD, "mesh/sphere_facets.xdmf", "r");
        // read_meshtags требует строго 4 аргумента:
        // 1. Сетка (*mesh)
        // 2. Имя Grid ("Grid")
        // 3. Имя Attribute (явно указываем std::string, чтобы компилятор понял тип)
        // 4. Базовый путь XPath ("/Xdmf/Domain")
        auto facet_tags = file_facets.read_meshtags(*mesh, "Grid", std::string("Grid"), "/Xdmf/Domain");
        file_facets.close();


        // Пространства функций
        // пространство для скорости
        auto basix_V = basix::create_element<U>(
            basix::element::family::P, basix::cell::type::tetrahedron, 2,
            basix::element::lagrange_variant::unset, basix::element::dpc_variant::unset, false);
        auto element_V = std::make_shared<fem::FiniteElement<U>>(basix_V, std::vector<std::size_t>{3}, false);
        auto V = std::make_shared<fem::FunctionSpace<U>>(fem::create_functionspace<U>(mesh, element_V));
        // пространство для давления
        auto basix_Q = basix::create_element<U>(
            basix::element::family::P, basix::cell::type::tetrahedron, 1,
            basix::element::lagrange_variant::unset, basix::element::dpc_variant::unset, false);
        auto element_Q = std::make_shared<fem::FiniteElement<U>>(basix_Q, std::vector<std::size_t>{}, false);
        auto Q = std::make_shared<fem::FunctionSpace<U>>(fem::create_functionspace<U>(mesh, element_Q));
        // задаем функции нужные для шагов
        auto u_n = std::make_shared<fem::Function<T>>(V);
        auto u_s = std::make_shared<fem::Function<T>>(V);
        auto u_next = std::make_shared<fem::Function<T>>(V);
        auto p_n = std::make_shared<fem::Function<T>>(Q);
        auto p_s = std::make_shared<fem::Function<T>>(Q);
        // задаем константы
        auto dt = std::make_shared<fem::Constant<T>>(0.0001);
        auto rho = std::make_shared<fem::Constant<T>>(1.0); // Плотность воздуха, кг/м^3
        auto mu = std::make_shared<fem::Constant<T>>(0.05);
        // костыль - мы создаем доп пространство для сохранения скоростей для визуализации, тк паравиев с параболами дружить отказался 
        auto basix_V_vis = basix::create_element<U>(
            basix::element::family::P, basix::cell::type::tetrahedron, 1,
            basix::element::lagrange_variant::unset, basix::element::dpc_variant::unset, false);
        auto element_V_vis = std::make_shared<fem::FiniteElement<U>>(basix_V_vis, std::vector<std::size_t>{3}, false);
        auto V_vis = std::make_shared<fem::FunctionSpace<U>>(fem::create_functionspace<U>(mesh, element_V_vis));
        
        // функция для сохранения скорости
        auto u_vis = std::make_shared<fem::Function<T>>(V_vis);
        u_vis->name = "Velocity"; 
        p_n->name = "Pressure"; 

        // граничные условия
        /*
        auto u_inlet = std::make_shared<fem::Function<T>>(V);
        u_inlet->interpolate([](auto x) -> std::pair<std::vector<T>, std::vector<std::size_t>> {
            std::size_t N = x.extent(1);
            std::vector<T> v(3 * N, 0.0);
            for (std::size_t i = 0; i < N; ++i) v[0 * N + i] = 1.0; 
            return {v, {3, N}};
        });
        */
        auto u_inlet = std::make_shared<fem::Function<T>>(V);
        u_inlet->interpolate([](auto x) -> std::pair<std::vector<T>, std::vector<std::size_t>> {
            std::size_t N = x.extent(1);
            std::vector<T> v(3 * N, 0.0);

            double W = 0.6; // ширина по Y
            double H = 0.6; // высота по Z
            double u_max = 1.5; // максимальная скорость в центре
            // парабола
            for (std::size_t i = 0; i < N; ++i){
                double y = x(1, i);
                double z = x(2, i);
                double y_norm = (y - 0.3) / (W / 2.0); 
                double z_norm = (z - 0.3) / (H / 2.0);
                double profile = u_max * (1.0 - y_norm * y_norm) * (1.0 - z_norm * z_norm);
                v[0 * N + i] = std::max(0.0, profile); 
            }
            return {v, {3, N}};
        });
        
        auto inlet_dofs = fem::locate_dofs_topological(*V->mesh()->topology_mutable(), *V->dofmap(), 2, facet_tags.find(1));
        fem::DirichletBC<T> bc_inlet(u_inlet, inlet_dofs);
        
        // на границах скорость зануляем
        auto u_zero = std::make_shared<fem::Function<T>>(V);
        u_zero->interpolate([](auto x) -> std::pair<std::vector<T>, std::vector<std::size_t>> {
            return {std::vector<T>(3 * x.extent(1), 0.0), {3, x.extent(1)}};
        });
        auto wall_facets = facet_tags.find(3);
        auto sphere_facets = facet_tags.find(4);
        std::vector<std::int32_t> noslip_facets;
        noslip_facets.insert(noslip_facets.end(), wall_facets.begin(), wall_facets.end());
        noslip_facets.insert(noslip_facets.end(), sphere_facets.begin(), sphere_facets.end());
        auto noslip_dofs = fem::locate_dofs_topological(*V->mesh()->topology_mutable(), *V->dofmap(), 2, noslip_facets);
        fem::DirichletBC<T> bc_noslip(u_zero, noslip_dofs);
        // устанавливаем нулевое давление на вылете из трубы
        auto p_zero = std::make_shared<fem::Function<T>>(Q);
        p_zero->interpolate([](auto x) -> std::pair<std::vector<T>, std::vector<std::size_t>> {
            return {std::vector<T>(x.extent(1), 0.0), {1, x.extent(1)}};
        });
        auto outlet_dofs = fem::locate_dofs_topological(*Q->mesh()->topology_mutable(), *Q->dofmap(), 2, facet_tags.find(2));
        fem::DirichletBC<T> bc_outlet(p_zero, outlet_dofs);

        // тест что сетка хорошо прочиталась
        std::cout << "Узлов на входе (Inlet): " << inlet_dofs.size() << std::endl;
        std::cout << "Узлов на стенках и сфере: " << noslip_dofs.size() << std::endl;

        // конектимся с питончиком
        auto mesh_const = std::shared_ptr<const mesh::Mesh<U>>(mesh);
        auto V_const = std::shared_ptr<const fem::FunctionSpace<U>>(V);
        auto Q_const = std::shared_ptr<const fem::FunctionSpace<U>>(Q);
        // создаем словари с нашими переменными
        std::map<std::string, std::shared_ptr<const fem::Function<T>>> coeffs = {
            {"u_n", std::shared_ptr<const fem::Function<T>>(u_n)},
            {"u_s", std::shared_ptr<const fem::Function<T>>(u_s)},
            {"p_n", std::shared_ptr<const fem::Function<T>>(p_n)},
            {"p_s", std::shared_ptr<const fem::Function<T>>(p_s)}
        };
        std::map<std::string, std::shared_ptr<const fem::Constant<T>>> consts = {
            {"dt", std::shared_ptr<const fem::Constant<T>>(dt)},
            {"rho", std::shared_ptr<const fem::Constant<T>>(rho)},
            {"mu", std::shared_ptr<const fem::Constant<T>>(mu)}
        };
        // вот тут делаем указатели на предсозданый UFL код
        fem::Form<T> a1 = fem::create_form<T>(*form_navier_stokes_chorin_a1, {V_const, V_const}, coeffs, consts, {}, {}, mesh_const);
        fem::Form<T> L1 = fem::create_form<T>(*form_navier_stokes_chorin_L1, {V_const}, coeffs, consts, {}, {}, mesh_const);

        fem::Form<T> a2 = fem::create_form<T>(*form_navier_stokes_chorin_a2, {Q_const, Q_const}, coeffs, consts, {}, {}, mesh_const);
        fem::Form<T> L2 = fem::create_form<T>(*form_navier_stokes_chorin_L2, {Q_const}, coeffs, consts, {}, {}, mesh_const);

        fem::Form<T> a3 = fem::create_form<T>(*form_navier_stokes_chorin_a3, {V_const, V_const}, coeffs, consts, {}, {}, mesh_const);
        fem::Form<T> L3 = fem::create_form<T>(*form_navier_stokes_chorin_L3, {V_const}, coeffs, consts, {}, {}, mesh_const);

        // cборка матриц для решателя
        la::petsc::Matrix A1(fem::petsc::create_matrix(a1), false);
        MatZeroEntries(A1.mat());
        fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A1.mat(), ADD_VALUES), a1, {bc_inlet, bc_noslip});
        MatAssemblyBegin(A1.mat(), MAT_FLUSH_ASSEMBLY); MatAssemblyEnd(A1.mat(), MAT_FLUSH_ASSEMBLY);
        fem::set_diagonal<T>(la::petsc::Matrix::set_fn(A1.mat(), INSERT_VALUES), *V, {bc_inlet, bc_noslip});
        MatAssemblyBegin(A1.mat(), MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A1.mat(), MAT_FINAL_ASSEMBLY);

        la::petsc::Matrix A2(fem::petsc::create_matrix(a2), false);
        MatZeroEntries(A2.mat());
        fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A2.mat(), ADD_VALUES), a2, {bc_outlet});
        MatAssemblyBegin(A2.mat(), MAT_FLUSH_ASSEMBLY); MatAssemblyEnd(A2.mat(), MAT_FLUSH_ASSEMBLY);
        fem::set_diagonal<T>(la::petsc::Matrix::set_fn(A2.mat(), INSERT_VALUES), *Q, {bc_outlet});
        MatAssemblyBegin(A2.mat(), MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A2.mat(), MAT_FINAL_ASSEMBLY);

        la::petsc::Matrix A3(fem::petsc::create_matrix(a3), false);
        MatZeroEntries(A3.mat());
        fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A3.mat(), ADD_VALUES), a3, {bc_inlet, bc_noslip});
        MatAssemblyBegin(A3.mat(), MAT_FLUSH_ASSEMBLY); MatAssemblyEnd(A3.mat(), MAT_FLUSH_ASSEMBLY);
        fem::set_diagonal<T>(la::petsc::Matrix::set_fn(A3.mat(), INSERT_VALUES), *V, {bc_inlet, bc_noslip});
        MatAssemblyBegin(A3.mat(), MAT_FINAL_ASSEMBLY); MatAssemblyEnd(A3.mat(), MAT_FINAL_ASSEMBLY);

        std::vector<std::int32_t> sphere_facets_indices = facet_tags.find(4);

        // 2. Вычисляем домены интегрирования (пары: индекс_ячейки, локальный_индекс_грани)
        std::vector<std::int32_t> integration_entities = dolfinx::fem::compute_integration_domains(
            dolfinx::fem::IntegralType::exterior_facet,
            *mesh->topology(),
            sphere_facets_indices
        );

        // 3. Формируем map в том виде, в котором его жестко требует компилятор
        std::map<dolfinx::fem::IntegralType, std::vector<std::pair<std::int32_t, std::span<const std::int32_t>>>> subdomains = {
            {
                dolfinx::fem::IntegralType::exterior_facet, 
                { {4, std::span<const std::int32_t>(integration_entities)} }
            }
        };

        // 4. Создаем форму
        fem::Form<T> drag_form = fem::create_form<T>(
            *form_navier_stokes_chorin_drag_force, 
            {},           // Пусто, так как нет пространств тестовых/пробных функций
            coeffs,       // Передаем u_n и p_n
            consts,       // Передаем mu
            subdomains,   // Теперь тип совпадает идеально!
            {}, 
            mesh_const
        );
        // настройки солеверов
        // матрица A1 несимметрична из-за конвекции
        // поэтому используем "bcgs" (BiCGSTAB — метод бисопряженных градиентов со стабилизацией)
        // предобуславливатель (штука, которая упрощает матрицу перед решением) — "bjacobi" (Блочный Якоби), он отлично работает параллельно - просто чото крутое вроде
        la::petsc::KrylovSolver solver1(MPI_COMM_WORLD); solver1.set_operator(A1.mat());
        la::petsc::options::set("ksp_type", "bcgs"); la::petsc::options::set("pc_type", "bjacobi"); solver1.set_from_options();
        // уравнение Пуассона дает красивую, идеально симметричную матрицу A2
        // поэтому здесь используется "cg" (Метод сопряженных градиентов — самый быстрый для таких матриц)
        // а вот предобуславливатель "hypre" — это алгебраический многосеточный метод (AMG) - шобы быстрее считалось
        la::petsc::KrylovSolver solver2(MPI_COMM_WORLD); solver2.set_operator(A2.mat());
        la::petsc::options::set("ksp_type", "bcgs"); la::petsc::options::set("pc_type", "hypre"); solver2.set_from_options();
        // матрица A3 тоже симметрична ("cg")
        // сама по себе эта система очень простая (в a3 у нас только inner(u,v)), поэтому тут достаточно самого базового предобуславливателя "jacobi"
        la::petsc::KrylovSolver solver3(MPI_COMM_WORLD); solver3.set_operator(A3.mat());
        la::petsc::options::set("ksp_type", "bcgs"); la::petsc::options::set("pc_type", "jacobi"); solver3.set_from_options();

        // подготавливаем вектора для решения СЛАУ
        la::Vector<T> b1(L1.function_spaces()[0]->dofmap()->index_map, L1.function_spaces()[0]->dofmap()->index_map_bs());
        la::Vector<T> b2(L2.function_spaces()[0]->dofmap()->index_map, L2.function_spaces()[0]->dofmap()->index_map_bs());
        la::Vector<T> b3(L3.function_spaces()[0]->dofmap()->index_map, L3.function_spaces()[0]->dofmap()->index_map_bs());
        // правильные типы с правильной памятью для решателя
        la::petsc::Vector _us(la::petsc::create_vector_wrap(*u_s->x()), false);
        la::petsc::Vector _ps(la::petsc::create_vector_wrap(*p_s->x()), false);
        la::petsc::Vector _unext(la::petsc::create_vector_wrap(*u_next->x()), false);
        // количество шагов
        int num_steps = 1000;
        std::cout << "Начинаем расчет CFD (Метод Чорина)..." << std::endl;

        // задаем имена функциям (это важно для ParaView)
        // u_n->name = "Velocity";
        // p_n->name = "Pressure";

        // создаем файл для записи
        io::XDMFFile file_out(MPI_COMM_WORLD, "results.xdmf", "w");
        file_out.write_mesh(*mesh);
        
        // записываем начальные условия (t = 0)
        u_vis->interpolate(*u_n); // интерполируем P2 -> P1 - делаем из парабол то, что схавает паравиев
        file_out.write_function(*u_vis, 0.0);
        file_out.write_function(*p_n, 0.0);

        for (int step = 1; step <= num_steps; ++step) {
            double t = step * dt->value[0];
            // пересобираем матрицу А1 (остальные можно не пересобирать они константы - они геометрические)
            MatZeroEntries(A1.mat());
            fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A1.mat(), ADD_VALUES), a1, {bc_inlet, bc_noslip});
            MatAssemblyBegin(A1.mat(), MAT_FLUSH_ASSEMBLY); 
            MatAssemblyEnd(A1.mat(), MAT_FLUSH_ASSEMBLY);
            fem::set_diagonal<T>(la::petsc::Matrix::set_fn(A1.mat(), INSERT_VALUES), *V, {bc_inlet, bc_noslip});
            MatAssemblyBegin(A1.mat(), MAT_FINAL_ASSEMBLY); 
            MatAssemblyEnd(A1.mat(), MAT_FINAL_ASSEMBLY);
            // ШАГ 1
            std::ranges::fill(b1.array(), 0.0); // очищаем старый вектор
            fem::assemble_vector(b1.array(), L1); // L1 вычисляем
            fem::apply_lifting(b1.array(), {a1}, {{bc_inlet, bc_noslip}}, {}, T(1)); // применение граничных условий, но как-то хитро - надо разобраться
            b1.scatter_rev(std::plus<T>()); // синхронизируем для многопоточности
            bc_inlet.set(b1.array(), std::nullopt);
            bc_noslip.set(b1.array(), std::nullopt);
            
            la::petsc::Vector _b1(la::petsc::create_vector_wrap(b1), false);
            solver1.solve(_us.vec(), _b1.vec()); // решаем
            u_s->x()->scatter_fwd(); // делимся с остальными решателями (ну когда многопоток)
            // дебаговый вывод, чтобы понять улетает ли солвер в 0 или бесконечность
            if (step % 10 == 1) {
                double norm_b1, norm_us;
                VecNorm(_b1.vec(), NORM_2, &norm_b1);  // Проверяем правую часть (должна быть > 0)
                VecNorm(_us.vec(), NORM_2, &norm_us);  // Проверяем решение (должно быть > 0)
    
                if (rank == 0) {
                std::cout << "--- ШАГ 1 ДЕБАГ ---" << std::endl;
                std::cout << "Норма вектора правой части (b1): " << norm_b1 << std::endl;
                std::cout << "Норма решения скорости (u_s): " << norm_us << std::endl;
                }
            }

            // ШАГ 2
            std::ranges::fill(b2.array(), 0.0);
            fem::assemble_vector(b2.array(), L2);
            fem::apply_lifting(b2.array(), {a2}, {{bc_outlet}}, {}, T(1));
            b2.scatter_rev(std::plus<T>());
            bc_outlet.set(b2.array(), std::nullopt);
            
            la::petsc::Vector _b2(la::petsc::create_vector_wrap(b2), false);
            solver2.solve(_ps.vec(), _b2.vec());
            p_s->x()->scatter_fwd();

            // ШАГ 3
            std::ranges::fill(b3.array(), 0.0);
            fem::assemble_vector(b3.array(), L3);
            fem::apply_lifting(b3.array(), {a3}, {{bc_inlet, bc_noslip}}, {}, T(1));
            b3.scatter_rev(std::plus<T>());
            bc_inlet.set(b3.array(), std::nullopt);
            bc_noslip.set(b3.array(), std::nullopt);
            
            la::petsc::Vector _b3(la::petsc::create_vector_wrap(b3), false);
            solver3.solve(_unext.vec(), _b3.vec());
            u_next->x()->scatter_fwd();

            // сохраняем результаты шага и переходим к следующему
            std::ranges::copy(u_next->x()->array(), u_n->x()->array().begin());
            std::ranges::copy(p_s->x()->array(), p_n->x()->array().begin());

            u_n->x()->scatter_fwd();
            p_n->x()->scatter_fwd();
            // интерполируем для записи и визуализации
            u_vis->interpolate(*u_n);
            file_out.write_function(*u_vis, t);
            file_out.write_function(*p_n, t);
            T local_drag = fem::assemble_scalar(drag_form);
            
            // мы считаем на нескольких ядрах, нужно сложить кусочки
            double local_drag_val = std::real(local_drag); // приводим к double на всякий случай, иначе говно какое-то
            double total_drag = 0.0;
            MPI_Reduce(&local_drag_val, &total_drag, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

            if (rank == 0) {
                std::cout << "Шаг времени: " << step << "/" << num_steps << " (t = " << t << ")" << std::endl;
                std::cout << "Лобовое сопротивление (Drag): " << total_drag << " Н" << std::endl;
                std::cout << "------------------------------------------------" << std::endl;
            }
        }
        file_out.close();
        std::cout << "Расчет успешно завершен!" << std::endl;
    }

    // возвращаем стандартный вывод обратно в консоль и закрываем файл
    if (rank == 0) {
        std::cout.rdbuf(original_cout_buf);
        log_file.close();
        std::cout << "Логи успешно сохранены в log.txt" << std::endl;
    }
    PetscFinalize();
    return 0;
}