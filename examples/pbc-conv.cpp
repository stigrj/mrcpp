#include "pbc.h"

int order = 5;
int min_scale = 0;
int max_depth = 5;
double prec = 1.0e-1;
double bond_length = 5.3;
std::string type = "ionic";

void apply_operator(int n,
                    const mrcpp::MultiResolutionAnalysis<D> &MRA,
                    mrcpp::FunctionTree<D> *f_tree,
                    std::vector<mrcpp::FunctionTree<D> *> &g_trees,
                    std::vector<double> &energies) {
    int r = 1;
    switch (-n-0) {
        case 0:
            r = 32;
            break;
        case 1:
            r = 16;
            break;
        case 2:
            r = 8;
            break;
        case 3:
            r = 4;
            break;
        case 4:
            r = 2;
            break;
        case 5:
            r = 1;
            break;
        case 6:
            r = 1;
            break;
        case 7:
            r = 1;
            break;
        case 8:
            r = 1;
            break;
        case 9:
            r = 1;
            break;
        case 10:
            r = 1;
            break;
    }
    //r = 1;

    mrcpp::Timer t0;
    mrcpp::PoissonOperator P(MRA, prec, n, r);
    t0.stop();
    auto *g_tree = new mrcpp::FunctionTree<D>(MRA);

    // Apply Poisson operator
    mrcpp::Timer t1;
    //mrcpp::Printer::setPrintLevel(10);
    //mrcpp::print::separator(0, ' ');
    //mrcpp::print::separator(0, ' ');
    //mrcpp::print::separator(0, ' ');
    //mrcpp::print::separator(0, '+');
    mrcpp::apply(prec, *g_tree, P, *f_tree);
    //mrcpp::Printer::setPrintLevel(0);
    t1.stop();

    // Compute properties
    auto g_int = g_tree->integrate();
    auto g_norm = std::sqrt(g_tree->getSquareNorm());
    auto d_norm = 0.0;
    auto d_energy = 0.0;
    auto n_energy = mrcpp::dot(*g_tree, *f_tree);
    if (energies.size() > 0) d_energy = n_energy - energies.back();

    double x = bond_length / 2.0;
    mrcpp::Coord<D> O{-3.0 * x, -x, -x};
    mrcpp::Coord<D> A{6.0 * x, 0.0, 0.0};
    mrcpp::Plotter<D> plt;
    plt.setOrigin(O);
    plt.setRange(A);

    // Plot charge
    {
        mrcpp::FunctionTree<D> f_tmp(MRA);
        mrcpp::copy_grid(f_tmp, *f_tree);
        mrcpp::copy_func(f_tmp, *f_tree);
        mrcpp::refine_grid(f_tmp, 1);
        std::string f_name = "f_tree-" + std::to_string(-n);
        plt.linePlot({10000}, f_tmp, f_name);
    }

    // Plot potential
    {
        mrcpp::FunctionTree<D> g_tmp(MRA);
        mrcpp::copy_grid(g_tmp, *g_tree);
        mrcpp::copy_func(g_tmp, *g_tree);
        mrcpp::refine_grid(g_tmp, 1);
        std::string g_name = "g_tree-" + std::to_string(-n);
        plt.linePlot({10000}, g_tmp, g_name);
    }

    // Plot difference
    if (g_trees.size() > 0) {
        mrcpp::FunctionTree<D> d_tree(MRA);
        mrcpp::build_grid(d_tree, *g_tree);
        mrcpp::build_grid(d_tree, *g_trees.front());
        mrcpp::add(-1.0, d_tree, 1.0, *g_tree, -1.0, *g_trees.back());
        d_norm = d_tree.getSquareNorm();
        mrcpp::refine_grid(d_tree, 1);
        std::string d_name = "d_tree-" + std::to_string(-n);
        plt.linePlot({10000}, d_tree, d_name);
    }

    // Print properties
    println(0, std::setw(4) << n <<
               std::setw(4) << r <<
               std::setw(10) << std::setprecision(1) << t0.elapsed() <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << std::setprecision(12) << n_energy <<
               std::setw(10) << std::setprecision(1) << d_energy <<
               std::setw(21) << std::setprecision(12) << g_norm <<
               std::setw(13) << std::setprecision(1) << d_norm);

    energies.push_back(n_energy);
    g_trees.push_back(g_tree);
}

int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    if (argc > 1) type = argv[1];
    if (argc > 2) prec = std::atof(argv[2]);
    if (argc > 3) order = std::atoi(argv[3]);
    if (argc > 4) min_scale = -std::atoi(argv[4]);
    if (argc > 5) MSG_ABORT("Too many arguments!");

    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::Printer::setWidth(95);
    mrcpp::print::environment(0);

    println(0, "type  : " << type);
    println(0, "prec  : " << prec);
    println(0, "order : " << order);
    println(0, "scale : " << min_scale);
    mrcpp::print::separator(0, ' ', 1);

    auto corner = std::array<int, D>{-1, -1, -1};
    auto boxes = std::array<int, D>{2, 2, 2};
    auto sfac = std::array<double, 3>{bond_length, bond_length, bond_length};
    auto world = mrcpp::BoundingBox<3>(0, corner, boxes, sfac, true);
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);
    MRA.print();

    mrcpp::FunctionTree<D> *f_tree = nullptr;
    std::vector<mrcpp::FunctionTree<D> *> g_trees;
    std::vector<double> energies;

    mrcpp::print::header(0, "Projecting charge density");
    println(0, std::setw(4) << "f" <<
               std::setw(10) << "t_proj" <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << "sq_norm" <<
               std::setw(21) << "integral");
    mrcpp::print::separator(0, '-');
    if (type == "ionic") f_tree = project_ionic(MRA);
    if (type == "covalent") f_tree = project_covalent(MRA);
    if (type == "monopole") f_tree = project_monopole(MRA);
    if (type == "dipole") f_tree = project_dipole(MRA);
    if (type == "quadrupole") f_tree = project_quadrupole(MRA);
    mrcpp::print::footer(0, timer, 2);

    mrcpp::print::header(0, "Applying Poisson operator");
    println(0, std::setw(4) << "n" <<
               std::setw(4) << "r" <<
               std::setw(10) << "t_build" <<
               std::setw(10) << "t_apply" <<
               std::setw(21) << "energy" <<
               std::setw(10) << "update" <<
               std::setw(21) << "sq_norm" <<
               std::setw(13) << "diff_norm");
    mrcpp::print::separator(0, '-');
    for (int n = 0; n >= min_scale; n--) apply_operator(n, MRA, f_tree, g_trees, energies);
    mrcpp::print::footer(0, timer, 2);

    delete f_tree;
    for (int i = 0; i < g_trees.size(); i++) delete g_trees[i];

    return 0;
}
