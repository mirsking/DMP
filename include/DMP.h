#pragma once

#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct sDMP{
	MatrixXd x, xd, xdd;//演示数据 
	VectorXd t_dem, t_run, x0, g, gm;
	double vm;
	unsigned int n;//时间向量的长度
	//初始化数据
	double D;
	double K;
	double tau,dt;
	double alpha;
	unsigned int n_w;
	VectorXd wc, wh;//权重函数时间度、中心和权重函数高度
	VectorXd s;
	VectorXd f;
	MatrixXd psi;
	MatrixXd w0,w1,w2;
};

class CDMP
{
public:
	CDMP(void);
	CDMP(double t_dem, double dt);
	~CDMP(void);
	int load(string DemPath);
	int diff();
	MatrixXd diff(MatrixXd x);
	void learnDMP();
	void Save2TXT(string path, MatrixXd x);
	sDMP dmpPara;

};