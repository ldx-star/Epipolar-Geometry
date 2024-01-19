//
// Created by ldx on 24-1-19.
//
#include<iostream>
#include<fstream>
#include<sstream>
#include<cassert>
#include<string>
#include<vector>
#include<opencv2/opencv.hpp>
#include<set>
#include<cstdlib>
#include<Eigen/Core>
#include<Eigen/SVD>
/***
 * 计算RANSAC采样成功所需的迭代次数  N = log(1-z) / log(1 - p^k)
 * @param p 内点率
 * @param k 拟合模型需要的样本个数，求解基础矩阵 k=8
 * @param z 预期采样成功概率
 * @return RANSAC采样成功所需的迭代次数
 */
int cal_ransac_iterations(double p, int k, double z = 0.99) {
    return static_cast<int>(round(std::log(1 - z) / std::log(1 - std::pow(p, k))));
}
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

std::vector<int> find_inters(const std::vector<std::pair<cv::Point2d,cv::Point2d>>& pair_points,const Eigen::Matrix<double,3,3>& F,const double inter_thresh){

}

int main() {
    srand(time(0));

    std::vector<std::pair<cv::Point2d, cv::Point2d>> pair_points;

    std::ifstream in("../correspondences.txt");
    assert(in.is_open());


    std::string line, word;
    int n_line = 0;
    while (getline(in, line)) {
        std::stringstream stream(line);
        if (n_line == 0) {
            // 第一行是匹配点对儿的数量
            int n_corrs = 0;
            stream >> n_corrs;
            pair_points.resize(n_corrs);
            n_line++;
        } else {
            double x1, y1, x2, y2;
            stream >> x1 >> y1 >> x2 >> y2;
            pair_points[n_line - 1].first.x = x1;
            pair_points[n_line - 1].first.y = y1;
            pair_points[n_line - 1].second.x = x2;
            pair_points[n_line - 1].second.y = y2;
        }
        n_line++;
    }

    //计算采样次数
    const float inter_ratio = 0.5; //内点率
    const int n_samples = 8;
    int n_iterations = cal_ransac_iterations(inter_ratio,n_samples);
    const double inter_thresh = 0.0015;//内点门线设置

    std::vector<int> best_inter; //最终估计内点
    std::cout << "RANSAC-F: Running for " << n_iterations
              << " iterations, threshold " << inter_thresh
              << "..." << std::endl;

    for(int i = 0; i < n_iterations; i++){
        // 随机找到8对儿不重复的点
        std::set<int> indexes;
        while(indexes.size() < 8){
            indexes.insert(rand() % pair_points.size());
        }
        std::set<int>::const_iterator iter = indexes.cbegin();
        Eigen::Matrix<double,3,8> pset1,pset2;
        for(int j = 0; j < 8 ; j++) {
            auto match = pair_points[*iter];
            pset1(0, j) = match.first.x;
            pset1(1, j) = match.first.y;
            pset1(2, j) = 1.0;

            pset2(0, j) = match.second.x;
            pset2(1, j) = match.second.y;
            pset2(2, j) = 1.0;
        }

        // 8点法估计基础矩阵
        Eigen::Matrix<double,3,3> F = fundamental_8_point(pset1,pset2);

        // 统计内点个数
        std::vector<int> inter_indexes = find_inters(pair_points,F,inter_thresh);




    }

}