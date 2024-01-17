//
// Created by ldx on 24-1-16.
//
#include<iostream>
#include<Eigen/Core>

class Camera {

public:
    Camera() {
        //归一化坐标，不考虑图像尺寸
        _c[0] = 0;
        _c[1] = 0;
    }

    Eigen::Vector2d projection(Eigen::Vector3d p3d) {
        /// 世界坐标 -> 相机坐标
        double xc = p3d[0] * _R[0] + p3d[1] * _R[1] + p3d[2] * _R[2] + _t[0];
        double yc = p3d[0] * _R[3] + p3d[1] * _R[4] + p3d[2] * _R[5] + _t[1];
        double zc = p3d[0] * _R[6] + p3d[1] * _R[7] + p3d[2] * _R[8] + _t[2];

        /// 相机坐标系 -> 像平面
        double x = _f * (xc / zc);
        double y = _f * (yc / zc);

        /// 径向畸变
        double r2 = x * x + y * y;
        double distort_ratio = 1 + _dist[0] * r2 + _dist[1] * r2 * r2;

        /// 像平面 -> 像素坐标
        Eigen::Vector2d p;
        p[0] = distort_ratio * x;
        p[1] = distort_ratio * y;
        return p;
    }

    Eigen::Vector3d pos_in_world() {
        /**
         *  相机中心在相机坐标系的坐标点是 （0，0，0）
         *  p_c = Rp_w + t
         *  p_w = -R^Tt (旋转矩阵 R^T = R^-1)
         */
        Eigen::Vector3d pos;
        pos[0] = _R[0] * _t[0] + _R[3] * _t[1] + _R[6] * _t[2];
        pos[1] = _R[1] * _t[0] + _R[4] * _t[1] + _R[7] * _t[2];
        pos[2] = _R[2] * _t[0] + _R[5] * _t[1] + _R[8] * _t[2];
        return pos;
    }

public:
    double _f;//焦距
    double _dist[2];//径向畸变参数
    double _c[2];//中心点坐标
    /**
     *      [ R[0] R[1] R[2] ]
     *      [ R[3] R[4] R[5] ]
     *      [ R[6] R[7] R[8] ]
     */
    double _R[9];//旋转矩阵
    double _t[3];//平移矩阵
};

int main() {
    Camera cam;
    cam._f = 0.920227;
    cam._dist[0] = -0.106599;
    cam._dist[1] = 0.104385;
    cam._t[0] = 0.0814358;
    cam._t[1] = 0.937498;
    cam._t[2] = -0.0887441;
    cam._R[0] = 0.999796;
    cam._R[1] = -0.0127375;
    cam._R[2] = 0.0156807;
    cam._R[3] = 0.0128557;
    cam._R[4] = 0.999894;
    cam._R[5] = -0.0073718;
    cam._R[6] = -0.0155846;
    cam._R[7] = 0.00757181;
    cam._R[8] = 0.999854;

    Eigen::Vector3d p3d;//世界坐标
    p3d << 1.36939, -1.17123, 7.04869;

    ///世界坐标的投影点
    Eigen::Vector2d p2d = cam.projection(p3d);
    std::cout << "projection coord:\n " << p2d[0] << " " << p2d[1] << std::endl;
    std::cout << "result should be:\n 0.208188 -0.035398\n\n";

    ///计算相机在世界坐标系的位置
    Eigen::Vector3d pos =cam.pos_in_world();
    std::cout<<"cam position in world is:\n "<< pos[0] << " " << pos[1] << " " << pos[2] <<std::endl;
    std::cout<<"result should be: \n -0.0948544 -0.935689 0.0943652\n\n";
    return 0;
}