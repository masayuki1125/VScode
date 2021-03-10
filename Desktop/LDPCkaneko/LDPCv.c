#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>
#include<complex.h>
#pragma warning(disable : 4996)
#include"MT.h"
#define Wc 2 //固定
#define Wr 4//2の倍数
#define N 1200
#define M (N*Wc/Wr)
#define Lmax 100
#define NOCULC -1

//i行目とj行目を入れ替える
void swap_row(int **H, int i,int j)
{
  int k,a;
  for(k=0;k<N;k++){
    a=H[i][k];
    H[i][k]=H[j][k];
    H[j][k]=a;
  }
}


//i列目とj列目を入れ替える
void swap_column(int **H,int i,int j)
{
  int k,a;
  for(k=0;k<M;k++){
    a=H[k][i];
    H[k][i]=H[k][j];
    H[k][j]=a;
  }
}

//H[i][x]列目を削除する
void delete_row(int** H, int i) {
	for (int j = i+1; j < M; j++) {
		for (int k = 0; k < N; k++) {
			H[j-1][k] = H[j][k];
		}
	}
	//最後の行を全部0にする
	for (int i = 0; i < N; i++) {
		H[M - 1][i] = 0;
	}
}

//i行目にj行目を足す
void add(int** H, int i, int j) {
	int k;
	for (k = 0; k < N; k++) {
		H[i][k] = H[i][k] ^ H[j][k];
	}
}

//転置行列の作成
void trancepose_matrix(int** A,int** B, int column, int row){
	for (int i = 0; i < column; i++) {
		for (int j = 0; j < row; j++) {
			B[j][i] = A[i][j];
		}
	}
}

void initialize_matrix(int **H,int column,int row)
{
	for (int i = 0; i < column; i++)
	{
		for (int j = 0; j < row; j++)
		{
			H[i][j] = (int)0;
		}
	}
}

void duplicate_matrix(int** A,int** B, int column, int row){
	for (int i = 0; i < column; i++){
		for(int j=0;j<row;j++){
			B[i][j] = A[i][j];
		}
	}
}

void print_matrix(int** C,int column,int row)
{
	for (int i1 = 0; i1 < column; i1++)
	{
		for (int i2 = 0; i2 < row; i2++)
		{
			printf("%d,",C[i1][i2]);
		}
		printf("\n");
	}
	printf("\n");
}

//gallagerの構成法
int make_matrix(int** H)
{   //main block
	for (int i = 0; i < M/Wc; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (j >= Wr * i && j <= (i + 1) * Wr - 1)
			{
					H[i][j] = 1;
			}
		}
	}

	for (int k = 0; k < (Wc-1); k++)
	{
		int tmp[N];
		//N個をランダムシャッフル
		for (int i = 0; i < N; i++)
		{
			tmp[i] = i;
		}

		for (int i = N; i > 1; i--)
		{
			int a = i - 1;
			int b = genrand_int32() % N;
			int tmp2 = tmp[a];
			tmp[a] = tmp[b];
			tmp[b] = tmp2;
		}

		//サブブロックに代入
		for (int j = 0; j < N; j++)
		{
			for (int i = (k+1)*M / Wc; i < (k+2)*M / Wc; i++)
			{
				H[i][j] = H[i - (k+1)*M / Wc][tmp[j]];
			}
		}
	}
	return 0;
}

//ガウス消去法
int gaussian_elimination(int **H,int** H1){
	//入れ替えた列を記憶しておく
	int A[M];
	for (int i = 0; i < M; i++) {
		A[i] = i;
	}

	for (int j = 0; j < M; j++) {
		//ノードを確定させる
		for (int i = M-j; i > 0; i--) {
			if (H[i - 1][N - 1-j] == 1) {
				swap_row(H, i - 1, M - 1-j);
				int tmp = A[i - 1];
				A[i - 1] = A[M - 1 - j];
				A[i - 1] = tmp;
				break;
			}
		}

		//一次従属な列をどうにかする
		if (H[M - 1 - j][N - 1 - j] != 1){

			for (int k = 0; k<N-1-j; k++) {
				if (H[M - 1 - j][k]==1) {
					swap_column(H, k, N - 1 - j);
					swap_column(H1, k, N - 1 - j);
					break;
				}
			}
		}

		if (H[M - 1 - j][N - 1 - j] != 1) {//それでもだめなら列の削除
			delete_row(H, M - 1 - j);
			delete_row(H1, A[M - 1 - j]);
		}
		
		//確定ノード以外の列の値を0にする
		for (int k = M - 1-j; k > 0; k--) {
			if (H[k - 1][N - 1-j] == 1) {
				if (k - 1 == M - 1 - j)
				{
					printf("error1!");
					exit(1);
				}
				add(H, k - 1, M - 1-j);
			}
		}
	}
	return 0;
}

//メッセージの入力と後退代入で符号化[M-1][N]
void encode(int** Ht, int** C) {
	for (int j = 0; j < (N-(M-1)); j++){
		for (int i = 0; i < (N - (M-1)); i++) {
			C[j][i] = genrand_int32() % 2;
		}

		for (int i = 0; i < (M-1); i++) { //i番目の数値を出したい。
			for (int k = 0; k < N-(M-1)+i; k++) {//N-M<=i<N
				C[j][i + N - (M-1)] ^= C[j][k] & Ht[k][i];
				
			}
		}
	}
}

void BPSK_AWGN(double Eb_N0db, int* c, double* z) {
	double Eb_N0 = pow(10.0, (double)(Eb_N0db) / 10);
	double sigma = sqrt(1 / (2.0 * Eb_N0));
	for (int i = 0; i < N; i++) {
		double u1 = genrand_real3(); //(0,1)乱数を作成
		double u2 = genrand_real3();
		double nI = sigma * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
		//double nQ = sigma * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
		if (c[i] == 0) {
			z[i] = -1 + nI;
		}
		else if(c[i]==1){
			z[i] = 1 + nI;
		}

		if (z[i] >= 0) {
			c[i] = 1;
		}
		else if (z[i] < 0) {
			c[i] = 0;
		}
	}
}

//BSCのビット誤り率
/*double probability(double Eb_N0db) {
	double Eb_N0 = pow(10.0, (double)(Eb_N0db) / 10);
	double sigma = sqrt(1 / (2.0 * Eb_N0));
	int ErrorMax = 1000;//誤り率を導出するときのサンプル点の数
	int Error = 0;
	int Count = 0;
	do{
		Count++;
		double u1 = genrand_real3(); //(0,1)乱数を作成
		double u2 = genrand_real3();
		double nI = sigma * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
		if (nI < -1)
			Error++;
	} while (Error < ErrorMax);
	return (double)Error / Count;
}*/

//sign関数
/*int sign(double x){
	if (x >= 0)
		return 1;
	else
		return -1;
}

//ギャラガーのf関数
double f(double x) {

	double a;
	if (x > 12.5) x = 12.5;
	else if (x < 1E-6) x = 1E-6;
	else if (x == 0) x = 1E-6;
	a = log((exp(x) + 1) / (exp(x) - 1));
	return a;

}*/

//対数尤度比
/*double lambda(int* c,int n,double p) {
	if (c[n] == 0) {
		return log((1 - p) / p);
	}
	else {
		return log(p / (1 - p));
	}
}*/

//C*Ht=Sの行列の積の計算
/*void decode(int* c,int** Ht,int* s) {
	for (int j = 0; j < (M-1); j++) {//row
		for (int k = 0; k < N; k++)
			s[j] ^= c[k] & Ht[k][j];
	}
}

//ベクトルSが0かどうか判別する
int judge(int* s) {
	for (int i = 0; i < M-1; i++) {
		if (s[i] != 0) {
			return 0;//false
		}
	}
	return 1;//true
}

//復号
void sum_product(int** H1, int** H1t , int* c,double Eb_N0db) {
	double Eb_N0 = pow(10.0, (double)(Eb_N0db) / 10);
	double sigma = sqrt(1 / (2.0 * Eb_N0));
	double** alpha;  //対数外部値比 [M-1][N]
	double** beta;   //対数事前値比 [M-1][N]
	double* temp1, *temp2, *temp3;
	int m,n,i;
	int l;      //反復回数
	int** A, ** B; //行列Hの1が立っている行インデックスと列インデックス  A[M-1][Wr] B[Wc][N]
	int* P;

	//動的メモリ確保

	alpha = malloc(sizeof(double*) * (M-1));
	beta = malloc(sizeof(double*) * (M-1));

	for (m = 0; m < (M-1); m++) {
		alpha[m] = malloc(sizeof(double) * N);
		beta[m] = malloc(sizeof(double) * N);
	}

	A = malloc(sizeof(int*) * (M-1));
	if (A==NULL) {
		exit(1);
	}
	for (m = 0; m < (M-1); m++) {
		A[m] = malloc(sizeof(int) * Wr);
	}

	B = malloc(sizeof(int*) * Wc);
	if (B == NULL) {
		exit(1);
	}
	for (m = 0; m < Wc; m++) {
		B[m] = malloc(sizeof(int) * N);
	}

	P = malloc(sizeof(int) * (M - 1));
	if (P == NULL) {
		exit(1);
	}
	for (int i = 0; i < M - 1; i++){
		P[i] = 0;
	}

	temp1 = malloc(sizeof(double) * (M-1));
	temp2 = malloc(sizeof(double) * (M-1));
	temp3 = malloc(sizeof(double) * N);


	//B(n)は正則でないので、値が入らないところを識別するために初期化をしなければならない
	for (int i = 0; i < Wc; i++){
		for (int j = 0; j < N; j++){
			B[i][j] = NOCULC;
		}
	}
	

	// A(m),B(n)を求める
	int detecA = 0;
	for (m = 0; m < M-1; m++) {
		i = 0;
		for (n = 0; n < N; n++) {
			if (H1[m][n] == 1) {
				A[m][i] = n;
				i++;
				detecA++;
			}
		}
	}

	int detecB = 0;
	for (n = 0; n < N; n++) {
		i = 0;
		for (m = 0; m < M-1; m++) {
			if (H1[m][n] == 1) {
				B[i][n] = m;
				i++;
				detecB++;
			}
		}
	}

	if (detecA != detecB || detecA != Wr * (M - 1)) {
		printf("error!377");
		exit(1);
	}

		//step1 初期化
	for (m = 0; m < (M-1); m++) {
		for (n = 0; n < N; n++) {
			beta[m][n] = 0;
			alpha[m][n] = 0;
		}
	}

	for (l = 0; l < Lmax; l++) {

		//対数尤度比
		double lambda[N];
		for (int i = 0; i < N; i++) {
			if (c[i] == 0) {
				lambda[i] = 2 / (sigma*sigma);
			}
			else {
				lambda[i] = -1 * 2 /(sigma*sigma);
			}
		}

		//step2 行処理

		//まず総和を求めておく
		double tmp1 = 0;//alphaの前項
		double tmp2 = 0;//alphaの後項

		for (m = 0; m < M-1; m++) {
			temp1[m] = 1;
			temp2[m] = 0;
			for (n = 0; n < Wr; n++) {
				temp1[m] *= sign(lambda[A[m][n]] + beta[m][A[m][n]]);
				temp2[m] += f(fabs(lambda[A[m][n]] + beta[m][A[m][n]]));
			}
		}

		//αを求める

		for (m = 0; m < M-1; m++) {
			for (n = 0; n < Wr; n++) {
				double a = lambda[A[m][n]] + beta[m][A[m][n]];
				alpha[m][A[m][n]] = temp1[m] / sign(a) * f(temp2[m] - f(fabs(a)));
			}

		}

		//step3 列処理

		//betaの総和をあらかじめ求めておく

		double tmp3 = 0;//alphaの総和
		for (n = 0; n < N; n++) {
			for (m = 0; m < Wc; m++) {
				if(B[m][n]!= NOCULC)
				beta[B[m][n]][n] += alpha[B[m][n]][n];
				
			}
		}

		//betaを求める
		for (n = 0; n < N; n++) {
			temp3[n] = 0;
			for (m = 0; m < Wc; m++) {
				if(B[m][n]!=NOCULC)
					temp3[n] += alpha[B[m][n]][n];
			}
		}

		for (n = 0; n < N; n++) {
			for (m = 0; m < Wc; m++) {
				if (B[m][n] != NOCULC)
				beta[B[m][n]][n] = temp3[n] - alpha[B[m][n]][n];
			}
		}
		//step4 一時推定語の計算

		double tmp = 0;
		for (n = 0; n < N; n++) {
			
			for (m = 0; m < Wc; m++) {
				if (B[m][n] != NOCULC)
				tmp += alpha[B[m][n]][n];
			}

			for (n = 0; n < N; n++) {
				if (sign(lambda[n]+tmp) == 1) c[n] = 0;
				else c[n] = 1;
			}
		}

		//step5 パリティ検査
		for (m = 0; m < M-1; m++) {
			P[m] = 0;
		}
		int e = 0;


		for (m = 0; m < M-1; m++) {
			for (n = 0; n < N; n++) {
				P[m] = P[m] ^ (c[n] * H1t[n][m]);
			}
			e = e | P[m];
		}

		if (e == 0) {
			printf("success!");
			break;
		}
		if (l = 99) {
			printf("error!");
		}
	}

	//メモリの解放
	for (m = 0; m < (M-1); m++) {
		free(alpha[m]);
		free(beta[m]);
	}
	free(alpha);
	free(beta);

	for (m = 0; m < (M-1); m++) {
		free(A[m]);
	}
	free(A);

	for (m = 0; m < Wc; m++) {
		free(B[m]);
	}
	free(B);
}*/

double f(double x) {

	double a;
	if (x > 12.5) x = 12.5;
	else if (x < 1E-6) x = 1E-6;
	else if (x == 0) x = 1E-6;
	a = log((exp(x) + 1) / (exp(x) - 1));
	return a;

}

//sign関数
double sign(double x) {

	if (x >= 0) return 1;
	if (x < 0)  return -1;


}




//sum-product復号法
void sum_product(int** H, double y[], int c[], double sigma) {


	double** alpha;  //対数外部値比 [M][N]
	double** beta;   //対数事前値比 [M][N]
	double a, b, d, asum;
	double* temp1, * temp2, * temp3;  //予め求めておいた総和、総積
	int m, n, k, i, j;  //ループ用変数
	int l;      //反復回数
	int e;
	int* P;  //パリティ検査に使う行列  [M]
	int** A, ** B; //行列Hの1が立っている行インデックスと列インデックス  A[M][Wr] B[Wc][N]
	int** Ht;  //行列Hの転置  [N][M]

	//動的メモリ確保

	alpha = malloc(sizeof(double*) * (M-1));
	beta = malloc(sizeof(double*) * (M-1));

	for (m = 0; m < M-1; m++) {
		alpha[m] = malloc(sizeof(double) * N);
		beta[m] = malloc(sizeof(double) * N);
	}


	Ht = malloc(sizeof(int*) * N);
	for (n = 0; n < N; n++) {
		Ht[n] = malloc(sizeof(int) * (M-1));
	}

	P = malloc(sizeof(int) * (M-1));

	A = malloc(sizeof(int*) * (M-1));
	for (m = 0; m < (M-1); m++) {
		A[m] = malloc(sizeof(int) * Wr);
	}

	B = malloc(sizeof(int*) * Wc);
	for (m = 0; m < Wc; m++) {
		B[m] = malloc(sizeof(int) * N);
	}

	temp1 = malloc(sizeof(double) * (M-1));
	temp2 = malloc(sizeof(double) * (M-1));
	temp3 = malloc(sizeof(double) * N);

	//Hを転置させる   

	for (n = 0; n < N; n++) {
		for (m = 0; m < (M-1); m++) {
			Ht[n][m] = H[m][n];
		}
	}

	//B(n)だけ初期化
	for (m = 0; m < Wc; m++) {
		for (n = 0; n < N; n++) {
			B[m][n] = NOCULC;
		}
	}

	// A(m),B(n)を求める
	for (m = 0; m < (M-1); m++) {
		i = 0;
		for (n = 0; n < N; n++) {
			if (H[m][n] == 1) {
				A[m][i] = n;
				i++;
			}
		}
	}

	for (n = 0; n < N; n++) {
		i = 0;
		for (m = 0; m < (M-1); m++) {
			if (H[m][n] == 1) {
				B[i][n] = m;
				i++;
			}
		}
	}


	//step1 初期化


	for (m = 0; m < (M-1); m++) {
		for (n = 0; n < N; n++) {
			beta[m][n] = 0;
			alpha[m][n] = 0;
		}
	}



	for (l = 0; l < Lmax; l++) {

		//A(m),B(n)を使って総和などを予め求めておく

		for (m = 0; m < (M-1); m++) {
			temp1[m] = 1;
			temp2[m] = 0;
			for (n = 0; n < Wr; n++) {
				temp1[m] *= sign(2 * y[A[m][n]] / (sigma * sigma) + beta[m][A[m][n]]);
				temp2[m] += f(fabs(2 * y[A[m][n]] / (sigma * sigma) + beta[m][A[m][n]]));

			}
		}



		//step2 行処理

		for (m = 0; m < (M-1); m++) {
			for (n = 0; n < Wr; n++) {
				a = 2 * y[A[m][n]] / (sigma * sigma) + beta[m][A[m][n]];
				alpha[m][A[m][n]] = temp1[m] / sign(a) * f(temp2[m] - f(fabs(a)));
			}

		}


		//αの総和を求めておく

		for (n = 0; n < N; n++) {
			temp3[n] = 0;
			for (m = 0; m < Wc; m++) {
				if (B[m][n] != NOCULC) {
					temp3[n] += alpha[B[m][n]][n];
				}
			}
		}


		//step3 列処理

		for (n = 0; n < N; n++) {
			for (m = 0; m < Wc; m++) {
				if (B[m][n] != NOCULC) {
					beta[B[m][n]][n] = temp3[n] - alpha[B[m][n]][n];
				}
			}
		}


		//step4 一時推定語の計算

		for (n = 0; n < N; n++) {
			if (sign(2 * y[n] / (sigma * sigma) + temp3[n]) == 1) c[n] = 1;
			else c[n] = 0;

		}


		//step5 パリティ検査


		for (m = 0; m < (M-1); m++) {
			P[m] = 0;
		}
		e = 0;


		for (m = 0; m < (M-1); m++) {
			for (n = 0; n < N; n++) {
				P[m] = P[m] ^ (c[n] * Ht[n][m]);
			}
			e = e | P[m];
		}

		if (e == 0) break;


	}


	/*     //行列の確認用 */
	/*     for(m=0;m<(M-2);m++){ */
	/*       for(n=0;n<N;n++){ */
	/* 	printf("%f",alpha[m][n]); */
	/*       } */
	/*       printf("\n"); */
	/*     } */
	/*     printf("\n"); */



	//メモリの解放
	for (n = 0; n < N; n++) {
		free(Ht[n]);
	}
	free(Ht);

	for (m = 0; m < (M-1); m++) {
		free(alpha[m]);
		free(beta[m]);
	}
	free(alpha);
	free(beta);
	free(P);

	for (m = 0; m < (M-1); m++) {
		free(A[m]);
	}
	free(A);

	for (m = 0; m < Wc; m++) {
		free(B[m]);
	}
	free(B);

	free(temp1);
	free(temp2);
	free(temp3);


}


int main()
{
	int** H;//[M][N]
	int** Ht;//[N][M-1]
	int** H1;//[M][N]
	int** H1t;//[N][M-1]
	int** C;//[N-M+1][N]
	int** C1;//[N-M+1][N]
	int ErrorMax = 1000;
	int CountAll = 0;//bit送信数
	int CountError = 0;
	double Eb_N0db;
	init_genrand((unsigned)time(NULL));

	FILE* fp;
	fp = fopen("test.txt", "w");
	if (fp == NULL)
	{
		printf("cannot open\n");
		exit(1);
	}

	for (int i = 0; i < 10; i++) {

	//動的メモリ確保
	H = malloc(sizeof(int*) * M);
	H1 = malloc(sizeof(int*) * M);
	if (H == NULL|| H1==NULL){
		exit(1);
	}
	for (int i = 0; i < M; i++) {
		H[i] = malloc(sizeof(int) * N);
		H1[i] = malloc(sizeof(int) * N);
	}

	Ht = malloc(sizeof(int*) * N);
	H1t = malloc(sizeof(int*) * N);
	if (Ht == NULL || H1t==NULL) {
		exit(1);
	}
	for (int i = 0; i < N; i++) {
		Ht[i] = malloc(sizeof(int) * (M-1));
		H1t[i] = malloc(sizeof(int) * (M-1));
	}

	C = malloc(sizeof(int*) * (N-M+1));
	C1 = malloc(sizeof(int*) * (N - M + 1));
	if (C == NULL || C1==NULL){
		exit(1);
	}
	for (int i = 0; i < (N-M+1); i++) {
		C[i] = malloc(sizeof(int) * N);
		C1[i] = malloc(sizeof(int) * N);
	}

	
		Eb_N0db = i ;
		double Eb_N0 = pow(10.0, (double)(Eb_N0db) / 10);
		double sigma = sqrt(1 / (2.0 * Eb_N0));
		CountAll = 0;
		CountError = 0;
		do {
			//行列の初期化
			initialize_matrix(H, M, N);
			initialize_matrix(Ht, N, M - 1);
			initialize_matrix(H1, M, N);
			initialize_matrix(H1t, N, M - 1);
			initialize_matrix(C, (N - M + 1), N);
			initialize_matrix(C1, (N - M + 1), N);

			

			//符号化
			make_matrix(H);
			duplicate_matrix(H, H1, M, N);
			gaussian_elimination(H, H1);
			trancepose_matrix(H, Ht, M - 1, N);
			trancepose_matrix(H1, H1t, M - 1, N);
			encode(Ht, C);
			duplicate_matrix(C, C1, N - M + 1, N);//エラー確認用


			//復号
			for (int m = 0; m < N-M + 1; m++) {
				
				int c[N];
				double z[N];

				for (int i = 0; i < N; i++) {
					c[i] = C[m][i];
				}

				BPSK_AWGN(Eb_N0db,c,z);

				sum_product(H1,z, c, sigma);

				for (int i = 0; i < N; i++) {
					if (C[m][i] != c[i])
						CountError++;
				}
				
			}
			CountAll +=N*(N - (M - 1));
		} while (CountError < ErrorMax);
		double BER = (double)CountError/CountAll;
		printf("%f,%lg\n", Eb_N0db, BER);
		fprintf(fp, "%f,%lg\n", Eb_N0db, BER);//x=Eb/N0[db], y=BERを出力

	//メモリの開放
	for (int i = 0; i < M; i++){
		free(H[i]);
		free(H1[i]);
	}
	free(H);
	free(H1);


	for (int i = 0; i < (N - M+1); i++){
		free(C[i]);
		free(C1[i]);
	}
	free(C);
	free(C1);

	for (int i = 0; i < N; i++) {
		free(Ht[i]);
		free(H1t[i]);
	}
	free(Ht);
	free(H1t);


	}
	fclose(fp);

	return 0;
}



