#pragma once

#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct sDMP{
	MatrixXd x, xd, xdd, t_dem, t_run;//演示数据
	MatrixXd y, yd, ydd, z, zd, zdd, x0, g;
	unsigned int n;//时间向量的长度
	double dt;
	//初始化数据
	double D;
	double K;
	double tau, g, A;
	double alpha;
	unsigned int n_w;
	MatrixXd t_nw, wc, wh;//权重函数时间度、中心和权重函数高度
	double s, psi, f;
};

class CDMP
{
public:
	CDMP(void);
	CDMP(double t_dem, double dt);
	~CDMP(void);
	int load(string DemPath);
	int Diff();
	void learnDMP();
	sDMP dmpPara;

};