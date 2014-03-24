#include "DMP.h"
#include <string.h>
#include "math.h"

CDMP::CDMP(void)
{
	
}

CDMP::CDMP(double t, double dt)
{
	dmpPara.dt = dt;
	dmpPara.t_dem.setLinSpaced(1.0 / dt + 1, 0, 1).adjoint();
	dmpPara.D = 15;
	dmpPara.K = pow(dmpPara.D,2)/4;
	dmpPara.tau = dmpPara.t_dem(dmpPara.t_dem.rows()-1) - dmpPara.t_dem(0);	
	dmpPara.alpha = 4;
	dmpPara.n_w = dmpPara.t_dem.size()/10;
	VectorXd t_w;
	t_w.setLinSpaced(dmpPara.n_w,0,1);
	dmpPara.wc = (-dmpPara.alpha / dmpPara.tau*t_w);
	dmpPara.wc=dmpPara.wc.array().exp();
	dmpPara.wh = ((diff(dmpPara.wc*0.65)).array().pow(2)).array().inverse();
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
		while (fin >> a >> b >> c)
			RowNum++;
		dmpPara.x.resize(RowNum, 3);
		fin.clear();
		fin.seekg(0);
		RowNum = 0;
		while (fin >> a >> b >> c)
		{
			dmpPara.x(RowNum, 0) = a;
			dmpPara.x(RowNum, 1) = b;
			dmpPara.x(RowNum, 2) = c;
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
int CDMP::diff()
{
	MatrixXd tmpMat = MatrixXd::Zero(dmpPara.x.rows(), dmpPara.x.cols());
	//x的差分
	tmpMat.middleRows(0, dmpPara.x.rows() - 1) = dmpPara.x.middleRows(1, dmpPara.x.rows() - 1) - dmpPara.x.middleRows(0, dmpPara.x.rows() - 1);;
	tmpMat.row(tmpMat.rows()-1) = tmpMat.row(tmpMat.rows()-2);
	dmpPara.xd = tmpMat / dmpPara.dt;

	tmpMat.middleRows(0, dmpPara.xd.rows() - 1) = dmpPara.xd.middleRows(1, dmpPara.xd.rows() - 1) - dmpPara.xd.middleRows(0, dmpPara.xd.rows() - 1);;
	tmpMat.row(tmpMat.rows() - 1) = tmpMat.row(tmpMat.rows() - 2);
	dmpPara.xdd = tmpMat/dmpPara.dt;
	return 0;
}

MatrixXd CDMP::diff(MatrixXd x)
{
	MatrixXd tmpMat = MatrixXd::Zero(x.rows(), x.cols());
	//x的差分
	tmpMat.middleRows(0, x.rows() - 1) = x.middleRows(1, x.rows() - 1) - x.middleRows(0, x.rows() - 1);;
	tmpMat.row(tmpMat.rows() - 1) = tmpMat.row(tmpMat.rows() - 2);
	return tmpMat;
}


void CDMP::learnDMP()
{
	
	int n_t = dmpPara.t_dem.size();
	dmpPara.s.setZero(n_t);
	dmpPara.s(0) = 1;
	for (int i = 1; i < n_t;i++)
		dmpPara.s(i) = dmpPara.s(i - 1) - dmpPara.alpha*dmpPara.s(i - 1)*dmpPara.dt / dmpPara.tau;
	//w(i)的学习
	dmpPara.x0 = dmpPara.x.row(0);
	MatrixXd M_tmp,X,Y,Z;
	dmpPara.psi.resize(dmpPara.wc.size(), n_t);
	X = -dmpPara.wh*M_tmp.setOnes(1, n_t);
	Y = M_tmp.setOnes(dmpPara.wc.size(), 1)*(dmpPara.s.adjoint());
	Z = dmpPara.wc*M_tmp.setOnes(1, n_t);
	dmpPara.psi = (X.array()*( (Y - Z).array().pow(2) ).array()).array().exp();//100*1001
	//dmpPara.psi = ((-dmpPara.wh*M_tmp.setOnes(1, n_t)).array()*(((M_tmp.setOnes(dmpPara.wc.size(), 1)*(dmpPara.s.adjoint()) - dmpPara.wc*M_tmp.setOnes(1, n_t)).array().pow(2)).array())).array().exp();
	for (int i = 0; i < dmpPara.x.cols(); i++)
	{
		dmpPara.vm = dmpPara.tau*dmpPara.xd(n_t-1,i);
		dmpPara.gm = dmpPara.x(n_t-1, i) - dmpPara.vm*dmpPara.tau + (dmpPara.vm*dmpPara.t_dem).array();
		X = (
			pow(dmpPara.tau, 2) * dmpPara.xdd.col(i)
			- dmpPara.K*(dmpPara.gm - dmpPara.x.col(i))
			- dmpPara.D*(dmpPara.vm*M_tmp.setOnes(n_t, 1) - dmpPara.tau*dmpPara.xd.col(i))
			);
		Y = (dmpPara.K*(dmpPara.gm.array() - dmpPara.x0(i))*(dmpPara.s.array()));
		dmpPara.f = (X+Y) / dmpPara.K;
		// locally weighted linear regression
		MatrixXd wDen, wNum;
		X = (M_tmp.setOnes(dmpPara.wc.size(), 1));
		Y = ((dmpPara.s.adjoint()).array().pow(2));
		wDen = dmpPara.psi.array()*(X*Y).array();
		M_tmp = wDen;
		wDen = M_tmp.rowwise().sum();
		Y = dmpPara.s.array()*dmpPara.f.array();
		wNum = dmpPara.psi.array()*(X*Y.adjoint()).array();
		M_tmp = wNum;
		wNum = M_tmp.rowwise().sum();
		switch (i)
		{
		case 0:
			dmpPara.w0.setZero(dmpPara.wc.size(), 1);
			dmpPara.w0 = wNum.array() / wDen.array();
			Save2TXT("w0.txt", dmpPara.w0);
			break;
		case 1:
			dmpPara.w1.setZero(dmpPara.wc.size(), 1);
			dmpPara.w1 = wNum.array() / wDen.array();
			Save2TXT("w1.txt", dmpPara.w1);
			break;
		case 2:
			dmpPara.w2.setZero(dmpPara.wc.size(), 1);
			dmpPara.w2= wNum.array() / wDen.array();
			Save2TXT("w2.txt", dmpPara.w2);
			break;
		default:
			break;
		}
	}
}


void CDMP::Save2TXT(string path, MatrixXd x)
{
	fstream fout(path, ios::out);
	if (fout)
	{
		for (int i = 0; i < x.rows(); i++)
		{
			for (int j = 0; j < x.cols(); j++)
			{
				if (j == x.cols() - 1)
					fout << x(i, j) << endl;
				else
					fout << x(i, j) << "\t";
			}
		}
	}
	else
		cerr << "无法打开文件" << endl;
}