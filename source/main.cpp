#include <iostream>
#include <fstream>
#include <Eigen\Dense>
#include "DMP.h"

using namespace std;
using namespace Eigen;


int main()
{
	CDMP dmp(1,0.001);
	//将演示数据存储到
	dmp.load("..\\..\\data\\dem.txt");
	dmp.diff();
	dmp.learnDMP();
	MatrixXd v = (dmp.dmpPara.xd.row(dmp.dmpPara.xd.rows() - 1)).adjoint();
	MatrixXd g = (dmp.dmpPara.x.row(dmp.dmpPara.x.rows() - 1)).adjoint();
	dmp.dmpPara.obstacle=(dmp.dmpPara.x.row(500)).adjoint();
	dmp.ReproductDMP(1, v, g);
	return 0;
}