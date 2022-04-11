#include "pbc.h"
#include "MRCPP/trees/TreeIterator.h"

int order = 5;
int min_scale = 0;
int max_depth = 25;
double prec = 1.0e-4;
double bond_length = 5.3;
std::string type = "ionic";

int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    if (argc > 1) type = argv[1];
    if (argc > 2) prec = std::atof(argv[2]);
    if (argc > 3) order = std::atoi(argv[3]);
    if (argc > 4) min_scale = std::atoi(argv[4]);
    if (argc > 5) MSG_ABORT("Too many arguments!");

    auto printlevel = 0;
    mrcpp::Printer::init(printlevel);
    mrcpp::Printer::setWidth(90);
    mrcpp::print::environment(0);

    println(0, "type  : " << type);
    println(0, "prec  : " << prec);
    println(0, "order : " << order);
    println(0, "scale : " << -min_scale);
    mrcpp::print::separator(0, ' ', 1);

    /*
    auto foo = setup_mra(-min_scale);
    mrcpp::PoissonOperator bar(foo, prec);
    return 0;
    {
    auto MRA = setup_mra(0);
    mrcpp::FunctionTree<D> foo(MRA);
    auto f = [] (const mrcpp::Coord<D> &r) { return 1.0; };
    mrcpp::project<D>(-1.0, foo, f);
    println(0, foo);
    { auto &parent = foo.getNode(mrcpp::NodeIndex<D>{-3, {-1,-1,-1}}); }
    { auto &parent = foo.getNode(mrcpp::NodeIndex<D>{-3, {-1,-1, 1}}); }
    { auto &parent = foo.getNode(mrcpp::NodeIndex<D>{-3, {-1, 1,-1}}); }
    { auto &parent = foo.getNode(mrcpp::NodeIndex<D>{-3, { 1,-1,-1}}); }
    { auto &parent = foo.getNode(mrcpp::NodeIndex<D>{-3, { 1, 1,-1}}); }
    { auto &parent = foo.getNode(mrcpp::NodeIndex<D>{-3, { 1,-1, 1}}); }
    { auto &parent = foo.getNode(mrcpp::NodeIndex<D>{-3, {-1, 1, 1}}); }
    { auto &parent = foo.getNode(mrcpp::NodeIndex<D>{-3, { 1, 1, 1}}); }
    println(0, foo);
    {
        Eigen::MatrixXd norms = Eigen::MatrixXd::Zero(20, 4);
        mrcpp::HilbertIterator<D> it(&foo);
        while (it.nextParent()) {
            auto &node = it.getNode();
            println(0, node);
            int n = node.getScale();
            double sfac = 1.0;//std::pow(2.0, n+1.0);
            double sq_norm = sfac*node.getSquareNorm();
            double s_norm = sfac*node.getScalingNorm();
            double w_norm = sfac*node.getWaveletNorm();
            norms(n+min_scale,0) = n;
            norms(n+min_scale,1) += sq_norm;
            norms(n+min_scale,2) += s_norm;
            norms(n+min_scale,3) += w_norm;
        }
        println(0, std::endl);
        println(0, norms);
    }
    }
    return 0;
    */

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
    if (type == "ionic") f_tree = project_ionic();
    if (type == "covalent") f_tree = project_covalent();
    if (type == "monopole") f_tree = project_monopole();
    if (type == "dipole") f_tree = project_dipole();
    if (type == "quadrupole") f_tree = project_quadrupole();
    mrcpp::print::footer(0, timer, 2);

    mrcpp::print::header(0, "Applying Poisson operator");

    // Setup environment
    auto MRA = setup_mra(-min_scale);
    mrcpp::PoissonOperator P(MRA, prec);
    auto *g_tree = new mrcpp::FunctionTree<D>(MRA);
    g_tree->setZero();
    mrcpp::refine_grid(*g_tree, 2);
    mrcpp::clear_grid(*g_tree);

    // Apply Poisson operator
    mrcpp::Timer t1;
    mrcpp::Printer::setPrintLevel(10);
    mrcpp::apply(-1.0, *g_tree, P, *f_tree);
    mrcpp::Printer::setPrintLevel(0);
    t1.stop();

    println(0, *g_tree);

    // Compute properties
    auto g_int = g_tree->integrate();
    auto g_norm = std::sqrt(g_tree->getSquareNorm());
    auto n_energy = mrcpp::dot(*g_tree, *f_tree);
    mrcpp::print::footer(0, timer, 2);

    delete f_tree;
    for (int i = 0; i < g_trees.size(); i++) delete g_trees[i];

    return 0;
}
