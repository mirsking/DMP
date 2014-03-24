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
	return 0;
}