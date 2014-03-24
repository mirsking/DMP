#pragma once

#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct sDMP{
	MatrixXd x, xd, xdd, t_dem, t_run;//��ʾ����
	MatrixXd y, yd, ydd, z, zd, zdd, x0, g;
	unsigned int n;//ʱ�������ĳ���
	double dt;
	//��ʼ������
	double D;
	double K;
	double tau, g, A;
	double alpha;
	unsigned int n_w;
	MatrixXd t_nw, wc, wh;//Ȩ�غ���ʱ��ȡ����ĺ�Ȩ�غ����߶�
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