#include <iostream>
#include <random>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>

using namespace std;
using namespace Eigen;

void calcKalmanFilter(double & xhat_new, double & P_new, double & g, double A, double b, double bu, double c, double rho_vv, double rho_ww, double u, double y, double xhat,double P)
{
	double	xhatm 		= A*xhat + bu * u;
	double	Pm 		= A * P * A.transpose() + rho_vv * b * b.transpose();
	g 			= (Pm * c)/(c.transpose() * P * c + rho_ww);
	xhat_new		= xhatm + g * (y - c.transpose() * xhatm);
	P_new			= (1 - g * c.transpose()) * Pm;
	return;
}

int main()
{
	//システム
	double		rho_vv 	= 0.01;	//システム雑音の分散
	double		u 	= 0.0;	//制御の値
	double		rho_ww 	= 10.0;	//観測雑音の分散
	const int	N 	= 1000;	//データ数
	double		delta_t = 0.05;	

	double		k_spring	= 1.0;	//ばね定数
	double		m_mass		= 1.0;	//質量
	std::ofstream file("normal_distribution.tsv");

	// 観測データ生成
	// 雑音信号の生成
	typedef double matrixType;
	typedef Matrix<matrixType, 2, N> Matrix_2xN_f;
	typedef Matrix<matrixType, 1, N> Matrix_1xN_f;
	typedef Matrix<matrixType, 2, 2> Matrix_2x2_f;
	typedef Matrix<matrixType, 2, 1> Matrix_2x1_f;
	
	Matrix_2x2_f	A;					//プラントシステム
	Matrix_2x2_f	P;	//共分散行列
	Matrix_2x1_f	b;	//システム雑音の行列
	Matrix_2x1_f	c;	//観測行列
	Matrix_2x1_f	bu;	//制御行列
	Matrix_2xN_f 	x;
	Matrix_2xN_f 	xhat;
	Matrix_1xN_f	y;
	Matrix_1xN_f	v;
	Matrix_1xN_f	w;
	Matrix_2x1_f	g;

	A << 	1, -delta_t,
		k_spring*delta_t/m_mass, 1;

	P <<	1.0,	0,
		0,	1.0;

	b <<	0,
	  	0;

      	c <<	1.0,
      		0;

	bu<< 	0,
		0;
	
	x.col(0) << 	10,
			0;
	//乱数生成
	std::random_device seed_gen;
	std::default_random_engine engine(seed_gen());
	
	// 平均0.0、標準偏差sqrt(rho_vv)で分布させる?
	std::normal_distribution<matrixType> dist_v(0.0, sqrt(rho_vv));
	for (std::size_t k = 0; k < N; ++k) {
		// 正規分布で乱数を生成する
		v(k) = dist_v(engine);
	}

	// 平均0.0、標準偏差sqrt(rho_ww)で分布させる?
	std::normal_distribution<matrixType> dist_w(0.0, sqrt(rho_ww));
	for (std::size_t k = 0; k < N; ++k) {
		// 正規分布で乱数を生成する
		w(k) = dist_w(engine);
	}
	
	
	//真値生成
	for (std::size_t k = 0; k < (N - 1); ++k) {
		x.col(k+1) = A * x.col(k) + b * v(k);
		y(k) = c.transpose() * x.col(k) + w(k);
	}
	

	for (std::size_t k = 0; k < (N - 1); ++k) {
		calcKalmanFilter(xhat.col(k+1), P, g, A, b, bu, c, rho_vv, rho_ww, u, y(k), xhat.col(k), P);
	}
	
	file << "time" << "," << "true x1" << "," << "true x2" << "," << "estimated xhat" << "," << "measured y" << "\t\n";
	for (std::size_t k = 0; k < (N - 1); ++k) {
		//cout << x.col(k).row(0) << x.col(k).row(1) << endl;
		file << (double)k*delta_t << "," << x.col(k).row(0) << "," << x.col(k).row(1) << "," << "," << y(k) << endl;
		//file << x.col(k) << "\t\n";
	}
//	for (std::size_t k = 0; k < (N - 1); ++k) {
//		//file << x.col(k+1) << "," << xhat.col(k+1) << "," << y(k+1) << "\t\n";
//		std::cout << x(1,k) << std::end;
//	}
}

