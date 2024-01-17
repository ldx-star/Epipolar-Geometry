//
// Created by ldx on 24-1-17.
//
#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>

Eigen::Matrix<double, 3, 3> fundamental_8_point(const Eigen::Matrix<double, 3, 8> &pset1, const Eigen::Matrix<double, 3, 8> &pset2) {
    Eigen::Matrix<double, 3, 3> F;
    Eigen::Matrix<double, 8, 9> A;
    for (int i = 0; i < 8; i++) {
        A(i, 0) = pset1(0, i) * pset2(0, i);
        A(i, 1) = pset1(1, i) * pset2(0, i);
        A(i, 2) = pset2(0, i);
        A(i, 3) = pset1(0, i) * pset2(1, i);
        A(i, 4) = pset1(1, i) * pset2(1, i);
        A(i, 5) = pset2(1, i);
        A(i, 6) = pset1(0, i);
        A(i, 7) = pset1(1, i);
        A(i, 8) = 1.0;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd1(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd vv = svd1.matrixV();

    Eigen::VectorXd f = vv.col(8);
    F << f[0], f[1], f[2], f[3], f[4], f[5], f[6], f[7],f[8];

    //奇异值约束
    Eigen::JacobiSVD<Eigen::MatrixXd> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::MatrixXd U,V;
    U = svd2.matrixU();
    V = svd2.matrixV();
    Eigen::VectorXd v_S = svd2.singularValues();
    Eigen::Matrix3d S =Eigen::Matrix3d::Zero();

    S(0,0) = v_S(0); S(1,1) = v_S(1);S(2,2) = 0;
    F = U * S * V.transpose();

    return F;
}

int main() {
    Eigen::Matrix<double, 3, 8> pset1;
    pset1(0, 0) = 0.180123;
    pset1(1, 0) = -0.156584;
    pset1(2, 0) = 1.0;
    pset1(0, 1) = 0.291429;
    pset1(1, 1) = 0.137662;
    pset1(2, 1) = 1.0;
    pset1(0, 2) = -0.170373;
    pset1(1, 2) = 0.0779329;
    pset1(2, 2) = 1.0;
    pset1(0, 3) = 0.235952;
    pset1(1, 3) = -0.164956;
    pset1(2, 3) = 1.0;
    pset1(0, 4) = 0.142122;
    pset1(1, 4) = -0.216048;
    pset1(2, 4) = 1.0;
    pset1(0, 5) = -0.463158;
    pset1(1, 5) = -0.132632;
    pset1(2, 5) = 1.0;
    pset1(0, 6) = 0.0801864;
    pset1(1, 6) = 0.0236417;
    pset1(2, 6) = 1.0;
    pset1(0, 7) = -0.179068;
    pset1(1, 7) = 0.0837119;
    pset1(2, 7) = 1.0;

    Eigen::Matrix<double, 3, 8> pset2;
    pset2(0, 0) = 0.208264;
    pset2(1, 0) = -0.035405;
    pset2(2, 0) = 1.0;
    pset2(0, 1) = 0.314848;
    pset2(1, 1) = 0.267849;
    pset2(2, 1) = 1.0;
    pset2(0, 2) = -0.144499;
    pset2(1, 2) = 0.190208;
    pset2(2, 2) = 1.0;
    pset2(0, 3) = 0.264461;
    pset2(1, 3) = -0.0404422;
    pset2(2, 3) = 1.0;
    pset2(0, 4) = 0.171033;
    pset2(1, 4) = -0.0961747;
    pset2(2, 4) = 1.0;
    pset2(0, 5) = -0.427861;
    pset2(1, 5) = 0.00896567;
    pset2(2, 5) = 1.0;
    pset2(0, 6) = 0.105406;
    pset2(1, 6) = 0.140966;
    pset2(2, 6) = 1.0;
    pset2(0, 7) = -0.15257;
    pset2(1, 7) = 0.19645;
    pset2(2, 7) = 1.0;

    Eigen::Matrix<double, 3, 3> fundamental;
    fundamental = fundamental_8_point(pset1, pset2);
    std::cout << "Fundamental matrix after singularity constraint is:\n " << fundamental << std::endl;

    std::cout << "Result should be: \n" << "-0.0315082 -0.63238 0.16121\n"
              << "0.653176 -0.0405703 0.21148\n"
              << "-0.248026 -0.194965 -0.0234573\n" << std::endl;
    return 0;
}