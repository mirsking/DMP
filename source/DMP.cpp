#include "DMP.h"

CDMP::CDMP(void)
{
	
}

CDMP::CDMP(double t, double dt)
{
	dmpPara.t_dem.setLinSpaced(1/dt+1,0,1);
	dmpPara.D = 15;
	dmpPara.K = pow(dmpPara.D,2)/4;
	dmpPara.tau = dmpPara.t_dem(dmpPara.t_dem.cols()-1) - dmpPara.t_dem(0);
	dmpPara.x0 = dmpPara.x;
}

CDMP::~CDMP(void)
{

}

//读入示教数据
int CDMP::load(string DemPath)
{
	fstream fin(DemPath);
	double a, b, c;
	int RowNum = 0;
	if (fin)
	{
		while (fin >> c >> b >> a)
			RowNum++;
		dmpPara.x.resize(RowNum, 1);
		dmpPara.y.resize(RowNum, 1);
		dmpPara.z.resize(RowNum, 1);
		fin.clear();
		fin.seekg(0);
		RowNum = 0;
		while (fin >> c >> b >> a)
		{
			dmpPara.x(RowNum, 0) = a;
			dmpPara.y(RowNum, 0) = b;
			dmpPara.x(RowNum, 0) = c;
			RowNum++;
		}
		fin.close();
	}
	else
	{
		cerr << "打开文件失败！" << endl;
		return -1;
	}
	return 0;
}


//对示教数据求差分
int CDMP::Diff()
{
	MatrixXd tmpMat = MatrixXd::Zero(dmpPara.x.rows(), dmpPara.x.cols());
	//x的差分
	tmpMat.middleRows(0, dmpPara.x.rows() - 1) = dmpPara.x.middleRows(1, dmpPara.x.rows() - 1) - dmpPara.x.middleRows(0, dmpPara.x.rows() - 1);;
	tmpMat.row(tmpMat.rows()-1) = tmpMat.row(tmpMat.rows()-2);
	dmpPara.xd = tmpMat;

	tmpMat.middleRows(0, dmpPara.xd.rows() - 1) = dmpPara.xd.middleRows(1, dmpPara.xd.rows() - 1) - dmpPara.xd.middleRows(0, dmpPara.xd.rows() - 1);;
	tmpMat.row(tmpMat.rows() - 1) = tmpMat.row(tmpMat.rows() - 2);
	dmpPara.xdd = tmpMat;
	
	//y的差分
	tmpMat.middleRows(0, dmpPara.y.rows() - 1) = dmpPara.y.middleRows(1, dmpPara.y.rows() - 1) - dmpPara.y.middleRows(0, dmpPara.y.rows() - 1);;
	tmpMat.row(tmpMat.rows() - 1) = tmpMat.row(tmpMat.rows() - 2);
	dmpPara.yd = tmpMat;

	tmpMat.middleRows(0, dmpPara.yd.rows() - 1) = dmpPara.yd.middleRows(1, dmpPara.yd.rows() - 1) - dmpPara.yd.middleRows(0, dmpPara.yd.rows() - 1);;
	tmpMat.row(tmpMat.rows() - 1) = tmpMat.row(tmpMat.rows() - 2);
	dmpPara.ydd = tmpMat;

	//z的差分
	tmpMat.middleRows(0, dmpPara.z.rows() - 1) = dmpPara.z.middleRows(1, dmpPara.z.rows() - 1) - dmpPara.z.middleRows(0, dmpPara.z.rows() - 1);;
	tmpMat.row(tmpMat.rows() - 1) = tmpMat.row(tmpMat.rows() - 2);
	dmpPara.zd = tmpMat;

	tmpMat.middleRows(0, dmpPara.zd.rows() - 1) = dmpPara.zd.middleRows(1, dmpPara.zd.rows() - 1) - dmpPara.zd.middleRows(0, dmpPara.zd.rows() - 1);;
	tmpMat.row(tmpMat.rows() - 1) = tmpMat.row(tmpMat.rows() - 2);
	dmpPara.zdd = tmpMat;

	return 0;
}


void CDMP::learnDMP()
{
	double dt = 0.001;

}
