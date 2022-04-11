#include "pbc.h"
#include "MRCPP/trees/MWTree.h"
#include "MRCPP/trees/MWNode.h"
#include "MRCPP/trees/TreeIterator.h"

int order = 5;
int min_scale = 0;
int max_depth = 30;
double prec = 1.0e-4;
double bond_length = 5.3;
std::string type = "ionic";

int main(int argc, char **argv) {
    auto timer = mrcpp::Timer();

    if (argc > 1) type = argv[1];
    if (argc > 2) prec = std::atof(argv[2]);
    if (argc > 3) order = std::atoi(argv[3]);
    if (argc > 4) min_scale = -std::atoi(argv[4]);
    if (argc > 5) MSG_ABORT("Too many arguments!");

    auto printlevel = 1;
    mrcpp::Printer::init(printlevel);
    mrcpp::Printer::setWidth(90);
    mrcpp::print::environment(0);

    println(0, "type  : " << type);
    println(0, "prec  : " << prec);
    println(0, "order : " << order);
    println(0, "scale : " << min_scale);
    mrcpp::print::separator(0, ' ');

    auto corner = std::array<int, D>{-1, -1, -1};
    auto boxes = std::array<int, D>{2, 2, 2};
    auto sfac = std::array<double, 3>{bond_length, bond_length, bond_length};
    auto world = mrcpp::BoundingBox<3>(0, corner, boxes, sfac, true);
    auto basis = mrcpp::InterpolatingBasis(order);
    auto MRA = mrcpp::MultiResolutionAnalysis<D>(world, basis, max_depth);
    MRA.print();

    mrcpp::print::header(0, "Projecting charge density");
    println(0, std::setw(4) << "f" <<
               std::setw(10) << "t_proj" <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << "sq_norm" <<
               std::setw(21) << "integral");
    mrcpp::print::separator(0, '-');
    mrcpp::FunctionTree<D> *f_tree = nullptr;
    if (type == "ionic") f_tree = project_ionic(MRA);
    if (type == "covalent") f_tree = project_covalent(MRA);
    if (type == "monopole") f_tree = project_monopole(MRA);
    if (type == "dipole") f_tree = project_dipole(MRA);
    if (type == "quadrupole") f_tree = project_quadrupole(MRA);
    mrcpp::print::footer(0, timer, 2);

    mrcpp::PoissonOperator P(MRA, prec, min_scale);
    auto *g_tree = new mrcpp::FunctionTree<D>(MRA);
    mrcpp::apply(prec, *g_tree, P, *f_tree);

    return 0;
}
