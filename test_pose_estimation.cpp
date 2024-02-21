//
// Created by ldx on 24-1-30.
//
#include<Eigen/Core>
#include<Eigen/SVD>
#include<Eigen/Dense>
#include<Eigen/SVD>
#include<iostream>
#include<vector>

/*第一个相机的内外参数*/
double f1 = 0.972222208;
/*第二个相机的内外参数*/
double f2 = 0.972222208;

//匹配点对
Eigen::Vector2d p1 = {0.18012331426143646, -0.15658402442932129};
Eigen::Vector2d p2 = {0.2082643061876297, -0.035404585301876068};

Eigen::Vector3d triangulation(const Eigen::Vector2d &p1, const Eigen::Vector2d &p2, const Eigen::Matrix3d &K1, const Eigen::Matrix3d &R1, const Eigen::Vector3d &t1, const Eigen::Matrix3d &K2,
                              const Eigen::Matrix3d &R2, const Eigen::Vector3d &t2) {
    //投影矩阵
    Eigen::Matrix<double,3,4> M1,M2;
    M1 << R1,t1; M2 << R2,t2;
    M1 = K1 * M1; M2 = K2 * M2;
    std::cout<<"M1: "<<M1<<std::endl;
    std::cout<<"M1 for fist pose should be\n"
             <<"0.972222 0 0 0\n"
             <<"0 0.972222 0 0\n"
             <<"0 0 1 0\n";

    std::cout<<"M2: "<<M2<<std::endl;
    std::cout<<"M2 for fist pose should be\n"
             <<" -0.957966 0.165734 -0.00707496 0.0774496\n"
             <<"0.164089 0.952816 0.102143 0.967341\n"
             <<"0.0250416 0.102292 -0.994439 0.0605768\n";

    //构造A矩阵
    Eigen::Matrix4d A;
    A(0,0) = p1(0)*M1(2,0)-M1(0,0);
    A(0,1) = p1(0)*M1(2,1)-M1(0,1);
    A(0,2) = p1(0)*M1(2,2)-M1(0,2);
    A(0,3) = p1(0)*M1(2,3)-M1(0,3);

    A(1,0) = p1(1)*M1(2,0)-M1(1,0);
    A(1,1) = p1(1)*M1(2,1)-M1(1,1);
    A(1,2) = p1(1)*M1(2,2)-M1(1,2);
    A(1,3) = p1(1)*M1(2,3)-M1(1,3);

    A(2,0) = p2(0)*M2(2,0)-M2(0,0);
    A(2,1) = p2(0)*M2(2,1)-M2(0,1);
    A(2,2) = p2(0)*M2(2,2)-M2(0,2);
    A(2,3) = p2(0)*M2(2,3)-M2(0,3);

    A(3,0) = p2(1)*M2(2,0)-M2(1,0);
    A(3,1) = p2(1)*M2(2,1)-M2(1,1);
    A(3,2) = p2(1)*M2(2,2)-M2(1,2);
    A(3,3) = p2(1)*M2(2,3)-M2(1,3);

    Eigen::Vector3d P;
    Eigen::JacobiSVD<Eigen::Matrix4d> svd(A,Eigen::ComputeFullU | Eigen::ComputeFullV);
    auto V = svd.matrixV();
    P(0)=V(0,3)/V(3,3);
    P(1)=V(1,3)/V(3,3);
    P(2)=V(2,3)/V(3,3);
    return P;
}

bool is_correct_pose(Eigen::Matrix3d const &R1, Eigen::Vector3d const &t1, Eigen::Matrix3d &R2, Eigen::Vector3d const &t2) {
    //内参矩阵
    Eigen::Matrix3d K1, K2;
    K1 = K2 = Eigen::Matrix3d::Zero();
    K1(0, 0) = K1(1, 1) = f1;
    K2(0, 0) = K2(1, 1) = f2;
    K1(2, 2) = K2(2, 2) = 1.0;

    Eigen::Vector3d p3p = triangulation(p1, p2, K1, R1, t1, K2, R2, t2);//通过三角化求三维坐标
    auto x1 = R1 * p3p +t1;
    auto x2 = R2 * p3p + t2;
    return x1(2) > 0.0f && x2(2) > 0.0f;
}

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

    std::vector<std::pair<Eigen::Matrix3d, Eigen::Vector3d>> poses(4);
    poses[0] = {R1, t1};
    poses[1] = {R1, t2};
    poses[0] = {R2, t1};
    poses[0] = {R2, t2};

    //第一个相机姿态可以设置为[I|0], 第二个相机的姿态[R|t]
    Eigen::Matrix3d R_1 = Eigen::Matrix3d::Identity();
    Eigen::Vector3d t_1 = Eigen::Vector3d::Zero();

    //判断位姿是否合理
    bool flags[4];
    for (int i = 0; i < 4; i++) {
        flags[i] = is_correct_pose(R_1, t_1, poses[i].first, poses[i].second);
    }
    if(flags[0] || flags[1] || flags[2] || flags[3]){
        for(int i = 0; i < 4; i++){
            if(flags[i]){
                R = poses[i].first;
                t = poses[i].second;
            }
        }
        return true;
    }
    return false;
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
        std::cout<<"Correct pose found!"<<std::endl;
        std::cout<<"R: "<<R<<std::endl;
        std::cout<<"t: "<<t<<std::endl;
    }
    std::cout<<"Result should be: \n";
    std::cout<<"R: \n"
             << "0.999827 -0.0119578 0.0142419\n"
             << "0.0122145 0.999762 -0.0180719\n"
             << "-0.0140224 0.0182427 0.999735\n";
    std::cout<<"t: \n"
             <<"0.0796625 0.99498 0.0605768\n";
    return 0;
}



