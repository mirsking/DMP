#pragma once

#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct sDMP{
	MatrixXd x, xd, xdd;//演示数据 
	VectorXd t_dem, t_run, x0, obstacle, g, gm, vm, v_point, vd_point,v_point_bk;
	//初始化数据
	double D;
	double K;
	double tau,dt;
	double alpha;
	unsigned int n_w;
	VectorXd wc, wh;//权重函数时间度、中心和权重函数高度
	VectorXd s;
	VectorXd f;
	double f_point;
	MatrixXd psi;
	MatrixXd w,w0,w1,w2;
};

class CDMP
{
public:
	CDMP(void);
	CDMP(double tdem, double dt);
	~CDMP(void);
	int load(string DemPath);
	int diff();
	MatrixXd diff(MatrixXd x);
	void learnDMP();
	void ReproductDMP(double trun, VectorXd v_m, VectorXd g_m);
	double psi_obstacle(VectorXd x,VectorXd v,int dim);
	void Save2TXT(string path, MatrixXd x);
	sDMP dmpPara;

};