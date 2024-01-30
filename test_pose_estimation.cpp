//
// Created by ldx on 24-1-30.
//
#include<Eigen/Core>
#include<Eigen/SVD>
#include<Eigen/Dense>
#include<iostream>
#include<vector>

/*第一个相机的内外参数*/
double f1 = 0.972222208;
/*第二个相机的内外参数*/
double f2 = 0.972222208;

bool calc_cam_pose(const Eigen::Matrix<double, 3, 3> &F, const Eigen::Matrix3d &K1, const Eigen::Matrix3d &K2, Eigen::Matrix3d &R, Eigen::Vector3d &t) {
    // 计算本质矩阵
    Eigen::Matrix3d E = K2.transpose() * F * K1;
    std::cout << "EssentialMatrix result is \n" << E << std::endl;
    std::cout << "EssentialMatrix should be: \n"
              << "-0.00490744 -0.0146139 0.34281\n"
              << "0.0212215 -0.000748851 -0.0271105\n"
              << "-0.342111 0.0315182 -0.00552454\n";
    Eigen::JacobiSVD<Eigen::Matrix3d> svd_E(E, Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto U = svd_E.matrixU();
    auto V = svd_E.matrixV();

//    std::cout << std::endl  << U << std::endl;
//    std::cout << std::endl<< V << std::endl;

    //定义W和Z
    Eigen::Matrix3d W = Eigen::Matrix3d::Zero();
    W(0, 1) = -1, W(1, 0) = 1, W(2, 2) = 1;
    Eigen::Matrix3d Z = Eigen::Matrix3d::Zero();
    Z(0, 1) = -1, Z(1, 0) = 1;

    Eigen::Matrix3d UWtVt = U * W.transpose() * V.transpose();
    Eigen::Matrix3d UWVt = U * W * V.transpose();

    Eigen::Matrix3d R1 = UWVt.determinant() * UWVt;
    Eigen::Matrix3d R2 = UWtVt.determinant() * UWtVt;
    Eigen::Vector3d t1 = U.col(2);
    Eigen::Vector3d t2 = -U.col(2);

    std::vector<std::pair<Eigen::Matrix3d,Eigen::Vector3d>> poses(4);
    poses[0] = {R1,t1};
    poses[1] = {R1,t2};
    poses[0] = {R2,t1};
    poses[0] = {R2,t2};

    //第一个相机姿态可以设置为[I|0], 第二个相机的姿态[R|t]
    

}

int main() {
    Eigen::Matrix<double, 3, 3> F;
    F << -0.0051918668202215884, -0.015460923969578466, 0.35260470328319654,
            0.022451443619913483, -0.00079225386526248181, -0.027885130552744289,
            -0.35188558059920161, 0.032418724757766811, -0.005524537443406155;
    Eigen::Matrix3d R;
    Eigen::Vector3d t;
    Eigen::Matrix3d K1 = Eigen::Matrix3d::Zero();
    Eigen::Matrix3d K2 = Eigen::Matrix3d::Zero();
    K1(0, 0) = K1(1, 1) = f1;
    K1(2, 2) = 1.0;
    K2(0, 0) = K2(1, 1) = f1;
    K2(2, 2) = 1.0;
//    std::cout << K1 << std::endl;
    if (calc_cam_pose(F, K1, K2, R, t)) {

    }

    return 0;
}



