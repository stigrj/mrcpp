#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Gaussians"
#include "MRCPP/Printer"
#include "MRCPP/Plotter"
#include "MRCPP/Timer"
#include "MRCPP/utils/math_utils.h"

constexpr int D = 3;

extern int order;
extern int min_scale;
extern int max_depth;
extern double prec;
extern double bond_length;

mrcpp::GaussExp<D> setup_monopole() {
    double x = bond_length / 2.0;
    mrcpp::GaussExp<D> f_exp;
    {
        // Setting up analytic Gaussian
        auto beta = 1.2;
        auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({ 0.0, 0.0, 0.0});
        f_exp.append(f_func);
    }
    return f_exp;
}

mrcpp::GaussExp<D> setup_dipole() {
    double x = bond_length / 2.0;
    double shift = bond_length / 4.0;
    mrcpp::GaussExp<D> f_exp;
    {
        // Setting up analytic Gaussian
        auto beta = 1.2;
        auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({ x, x, x + shift});
        f_exp.append(f_func);
    }
    {
        // Setting up analytic Gaussian
        auto beta = 0.6;
        auto alpha = -std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({ x, x, x - shift});
        f_exp.append(f_func);
    }
    return f_exp;
}

mrcpp::GaussExp<D> setup_quadrupole() {
    double x = bond_length / 2.0;
    double shift = 0.0;
    mrcpp::GaussExp<D> f_exp;
    {
        // Setting up analytic Gaussian
        auto beta = 1.2;
        auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({ 0.0, x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({ 0.0,-x + shift,-x + shift});
        f_exp.append(f_func);
    }
    {
        // Setting up analytic Gaussian
        auto beta = 0.6;
        auto alpha = -std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({ 0.0,-x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({ 0.0, x + shift,-x + shift});
        f_exp.append(f_func);
    }
    return f_exp;
}

mrcpp::GaussExp<D> setup_ionic() {
    double x = bond_length / 2.0;
    double shift = x / 2.0;
    mrcpp::GaussExp<D> f_exp;
    {
        // Setting up analytic Gaussian
        auto beta = 1.2;
        auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({ x + shift, x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({ x - shift, x - shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({ x - shift, x + shift, x - shift});
        f_exp.append(f_func);
        f_func.setPos({ x + shift, x - shift, x - shift});
        f_exp.append(f_func);
    }
    {
        // Setting up analytic Gaussian
        auto beta = 0.6;
        auto alpha = -std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({ x - shift, x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({ x + shift, x - shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({ x + shift, x + shift, x - shift});
        f_exp.append(f_func);
        f_func.setPos({ x - shift, x - shift, x - shift});
        f_exp.append(f_func);
    }
    return f_exp;
}

mrcpp::GaussExp<D> setup_covalent() {
    double x = bond_length / 2.0;
    double shift = 0.0;
    mrcpp::GaussExp<D> f_exp;
    {
        // Setting up analytic Gaussian
        auto beta = 1.2;
        auto alpha = std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({ x + shift, x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({-x + shift, x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({ x + shift,-x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({ x + shift, x + shift,-x + shift});
        f_exp.append(f_func);
        f_func.setPos({-x + shift,-x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({-x + shift, x + shift,-x + shift});
        f_exp.append(f_func);
        f_func.setPos({ x + shift,-x + shift,-x + shift});
        f_exp.append(f_func);
        f_func.setPos({-x + shift,-x + shift,-x + shift});
        f_exp.append(f_func);
    }
    {
        // Setting up analytic Gaussian
        auto beta = 0.6;
        auto alpha = -std::pow(beta / mrcpp::pi, 3.0 / 2.0);
        auto f_func = mrcpp::GaussFunc<D>(beta, alpha);
        f_func.setPos({ x + shift, x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({-x + shift, x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({ x + shift,-x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({ x + shift, x + shift,-x + shift});
        f_exp.append(f_func);
        f_func.setPos({-x + shift,-x + shift, x + shift});
        f_exp.append(f_func);
        f_func.setPos({-x + shift, x + shift,-x + shift});
        f_exp.append(f_func);
        f_func.setPos({ x + shift,-x + shift,-x + shift});
        f_exp.append(f_func);
        f_func.setPos({-x + shift,-x + shift,-x + shift});
        f_exp.append(f_func);
    }
    return f_exp;
}

mrcpp::FunctionTree<D> * project_monopole(const mrcpp::MultiResolutionAnalysis<D> &MRA) {
    // Setup environment
    auto *f_tree = new mrcpp::FunctionTree<D>(MRA);

    // Setup density
    auto s_fac = std::array<double, 3>{bond_length, bond_length, bond_length};
    auto f_exp = setup_monopole();
    auto p_exp = f_exp.periodify(s_fac, 8.0);

    // Project function
    mrcpp::Timer t1;
    mrcpp::build_grid(*f_tree, p_exp);
    mrcpp::project(prec, *f_tree, p_exp);
    t1.stop();

    // Compute properties
    auto f_int = f_tree->integrate();
    auto f_norm = std::sqrt(f_tree->getSquareNorm());

    // Print properties
    println(0, std::setw(4) << "mono" <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << std::setprecision(12) << f_norm <<
               std::setw(21) << std::setprecision(12) << f_int);
    return f_tree;
}

mrcpp::FunctionTree<D> * project_dipole(const mrcpp::MultiResolutionAnalysis<D> &MRA) {
    // Setup environment
    auto *f_tree = new mrcpp::FunctionTree<D>(MRA);

    // Setup density
    auto s_fac = std::array<double, 3>{1.0 * bond_length, 1.0 * bond_length, 1.0 * bond_length};
    auto f_exp = setup_dipole();
    println (0, f_exp);
    auto p_exp = f_exp.periodify(s_fac, 8.0);

    // Project function
    mrcpp::Timer t1;
    mrcpp::build_grid(*f_tree, p_exp);
    mrcpp::project(prec, *f_tree, p_exp);
    t1.stop();

    // Compute properties
    auto f_int = f_tree->integrate();
    auto f_norm = std::sqrt(f_tree->getSquareNorm());

    // Print properties
    println(0, std::setw(4) << "di" <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << std::setprecision(12) << f_norm <<
               std::setw(21) << std::setprecision(12) << f_int);
    return f_tree;
}

mrcpp::FunctionTree<D> * project_quadrupole(const mrcpp::MultiResolutionAnalysis<D> &MRA) {
    // Setup environment
    auto *f_tree = new mrcpp::FunctionTree<D>(MRA);

    // Setup density
    auto s_fac = std::array<double, 3>{2.0 * bond_length, 2.0 * bond_length, 2.0 * bond_length};
    auto f_exp = setup_quadrupole();
    auto p_exp = f_exp.periodify(s_fac, 8.0);

    // Project function
    mrcpp::Timer t1;
    mrcpp::build_grid(*f_tree, p_exp);
    mrcpp::project(prec, *f_tree, p_exp);
    t1.stop();

    // Compute properties
    auto f_int = f_tree->integrate();
    auto f_norm = std::sqrt(f_tree->getSquareNorm());

    // Print properties
    println(0, std::setw(4) << "quad" <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << std::setprecision(12) << f_norm <<
               std::setw(21) << std::setprecision(12) << f_int);
    return f_tree;
}

mrcpp::FunctionTree<D> * project_ionic(const mrcpp::MultiResolutionAnalysis<D> &MRA) {
    // Setup environment
    auto *f_tree = new mrcpp::FunctionTree<D>(MRA);

    // Setup density
    auto s_fac = std::array<double, 3>{1.0 * bond_length, 1.0 * bond_length, 1.0 * bond_length};
    auto f_exp = setup_ionic();
    auto p_exp = f_exp.periodify(s_fac, 8.0);

    // Project function
    mrcpp::Timer t1;
    mrcpp::build_grid(*f_tree, p_exp);
    mrcpp::project(prec, *f_tree, p_exp);
    t1.stop();

    // Compute properties
    auto f_int = f_tree->integrate();
    auto f_norm = std::sqrt(f_tree->getSquareNorm());

    // Print properties
    println(0, std::setw(4) << "ion" <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << std::setprecision(12) << f_norm <<
               std::setw(21) << std::setprecision(12) << f_int);
    return f_tree;
}

mrcpp::FunctionTree<D> * project_covalent(const mrcpp::MultiResolutionAnalysis<D> &MRA) {
    // Setup environment
    auto *f_tree = new mrcpp::FunctionTree<D>(MRA);

    // Setup density
    auto s_fac = std::array<double, 3>{2.0 * bond_length, 2.0 * bond_length, 2.0 * bond_length};
    auto f_exp = setup_covalent();
    auto p_exp = f_exp.periodify(s_fac, 8.0);

    // Project function
    mrcpp::Timer t1;
    mrcpp::build_grid(*f_tree, p_exp);
    mrcpp::project(prec, *f_tree, p_exp);
    t1.stop();

    // Compute properties
    auto f_int = f_tree->integrate();
    auto f_norm = std::sqrt(f_tree->getSquareNorm());

    // Print properties
    println(0, std::setw(4) << "cov" <<
               std::setw(10) << std::setprecision(1) << t1.elapsed() <<
               std::setw(21) << " " <<
               std::setw(10) << " " <<
               std::setw(21) << std::setprecision(12) << f_norm <<
               std::setw(21) << std::setprecision(12) << f_int);
    return f_tree;
}
