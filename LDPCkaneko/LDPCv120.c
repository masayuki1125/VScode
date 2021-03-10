#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include<math.h>
#include<time.h>
#include<complex.h>
#pragma warning(disable : 4996)
#include"MT.h"
#define Wc 2 //�Œ�
#define Wr 4//2�̔{��
#define N 120
#define M (N*Wc/Wr)
#define Lmax 100
#define NOCULC -1
#define Mout (M-Wc+1)

//
void swap_row(int **H, int i,int j)
{
  int k,a;
  for(k=0;k<N;k++){
    a=H[i][k];
    H[i][k]=H[j][k];
    H[j][k]=a;
  }
}

//i��ڂ�j��ڂ����ւ���
void swap_column(int **H,int i,int j)
{
  int k,a;
  for(k=0;k<M;k++){
    a=H[k][i];
    H[k][i]=H[k][j];
    H[k][j]=a;
  }
}

//H[i][x]��ڂ��폜����
void delete_row(int** H, int i) {
	for (int j = i+1; j < M; j++) {
		for (int k = 0; k < N; k++) {
			H[j-1][k] = H[j][k];
		}
	}
	//�Ō�̍s��S��0�ɂ���
	for (int i = 0; i < N; i++) {
		H[M - 1][i] = 0;
	}
}

//i�s�ڂ�j�s�ڂ𑫂�
void add(int** H, int i, int j) {
	int k;
	for (k = 0; k < N; k++) {
		H[i][k] = H[i][k] ^ H[j][k];
	}
}

//�]�u�s��̍쐬
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

//gallager�̍\���@
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
		//N�������_���V���b�t��
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

		//�T�u�u���b�N�ɑ��
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

//�K�E�X�����@
int gaussian_elimination(int **H,int** H1){
	//����ւ�������L�����Ă���
	int A[M];
	for (int i = 0; i < M; i++) {
		A[i] = i;
	}

	for (int j = 0; j < M; j++) {
		//�m�[�h���m�肳����
		for (int i = M-j; i > 0; i--) {
			if (H[i - 1][N - 1-j] == 1) {
				swap_row(H, i - 1, M - 1-j);
				int tmp = A[i - 1];
				A[i - 1] = A[M - 1 - j];
				A[i - 1] = tmp;
				break;
			}
		}

		//�ꎟ�]���ȗ���ǂ��ɂ�����
		if (H[M - 1 - j][N - 1 - j] != 1){

			for (int k = 0; k<N-1-j; k++) {
				if (H[M - 1 - j][k]==1) {
					swap_column(H, k, N - 1 - j);
					swap_column(H1, k, N - 1 - j);
					break;
				}
			}
		}

		if (H[M - 1 - j][N - 1 - j] != 1) {//����ł����߂Ȃ��̍폜
			delete_row(H, M - 1 - j);
			delete_row(H1, A[M - 1 - j]);
		}
		
		//�m��m�[�h�ȊO�̗�̒l��0�ɂ���
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

//���b�Z�[�W�̓��͂ƌ�ޑ���ŕ�����[M-1][N]
void encode(int** Ht, int** C) {
	for (int j = 0; j < (N-(M-1)); j++){
		for (int i = 0; i < (N - (M-1)); i++) {
			C[j][i] = genrand_int32() % 2;
		}

		for (int i = 0; i < (M-1); i++) { //i�Ԗڂ̐��l���o�������B
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
		double u1 = genrand_real3(); //(0,1)�������쐬
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

double f(double x) {

	double a;
	if (x > 12.5) x = 12.5;
	else if (x < 1E-6) x = 1E-6;
	else if (x == 0) x = 1E-6;
	a = log((exp(x) + 1) / (exp(x) - 1));
	return a;

}

//sign�֐�
double sign(double x) {

	if (x >= 0) return 1;
	if (x < 0)  return -1;
	exit(1);
}

//sum-product�����@
void sum_product(int** H, double y[], int c[], double sigma) {


	double** alpha;  //�ΐ��O���l�� [M][N]
	double** beta;   //�ΐ����O�l�� [M][N]
	double a, b, d, asum;
	double* temp1, * temp2, * temp3;  //�\�ߋ��߂Ă��������a�A����
	int m, n, k, i, j;  //���[�v�p�ϐ�
	int l;      //������
	int e;
	int* P;  //�p���e�B�����Ɏg���s��  [M]
	int** A, ** B; //�s��H��1�������Ă���s�C���f�b�N�X�Ɨ�C���f�b�N�X  A[M][Wr] B[Wc][N]
	int** Ht;  //�s��H�̓]�u  [N][M]

	//���I�������m��

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

	//H��]�u������   

	for (n = 0; n < N; n++) {
		for (m = 0; m < (M-1); m++) {
			Ht[n][m] = H[m][n];
		}
	}

	//B(n)����������
	for (m = 0; m < Wc; m++) {
		for (n = 0; n < N; n++) {
			B[m][n] = NOCULC;
		}
	}

	// A(m),B(n)�����߂�
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


	//step1 ������
	for (m = 0; m < (M-1); m++) {
		for (n = 0; n < N; n++) {
			beta[m][n] = 0;
			alpha[m][n] = 0;
		}
	}



	for (l = 0; l < Lmax; l++) {

		//A(m),B(n)���g���đ��a�Ȃǂ�\�ߋ��߂Ă���

		for (m = 0; m < (M-1); m++) {
			temp1[m] = 1;
			temp2[m] = 0;
			for (n = 0; n < Wr; n++) {
				temp1[m] *= sign(2 * y[A[m][n]] / (sigma * sigma) + beta[m][A[m][n]]);
				temp2[m] += f(fabs(2 * y[A[m][n]] / (sigma * sigma) + beta[m][A[m][n]]));

			}
		}



		//step2 �s����
		for (m = 0; m < (M-1); m++) {
			for (n = 0; n < Wr; n++) {
				a = 2 * y[A[m][n]] / (sigma * sigma) + beta[m][A[m][n]];
				alpha[m][A[m][n]] = temp1[m] / sign(a) * f(temp2[m] - f(fabs(a)));
			}

		}


		//���̑��a�����߂Ă���

		for (n = 0; n < N; n++) {
			temp3[n] = 0;
			for (m = 0; m < Wc; m++) {
				if (B[m][n] != NOCULC) {
					temp3[n] += alpha[B[m][n]][n];
				}
			}
		}


		//step3 �񏈗�

		for (n = 0; n < N; n++) {
			for (m = 0; m < Wc; m++) {
				if (B[m][n] != NOCULC) {
					beta[B[m][n]][n] = temp3[n] - alpha[B[m][n]][n];
				}
			}
		}


		//step4 �ꎞ�����̌v�Z

		for (n = 0; n < N; n++) {
			if (sign(2 * y[n] / (sigma * sigma) + temp3[n]) == 1) c[n] = 1;
			else c[n] = 0;

		}


		//step5 �p���e�B����


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

	//�������̉��
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
	int CountAll = 0;//bit���M��
	int CountError = 0;
	double Eb_N0db;
	init_genrand((unsigned)time(NULL));

	FILE* fp;
	fp = fopen("120.txt", "w");
	if (fp == NULL)
	{
		printf("cannot open\n");
		exit(1);
	}

	//���I�������m��
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

	for (int i = 0; i < 26; i++) {
		Eb_N0db = i*0.2;
		double Eb_N0 = pow(10.0, (double)(Eb_N0db) / 10);
		double sigma = sqrt(1 / (2.0 * Eb_N0));
		CountAll = 0;
		CountError = 0;
		do {
			//�s��̏�����
			initialize_matrix(H, M, N);
			initialize_matrix(Ht, N, M - 1);
			initialize_matrix(H1, M, N);
			initialize_matrix(H1t, N, M - 1);
			initialize_matrix(C, (N - M + 1), N);
			initialize_matrix(C1, (N - M + 1), N);

			//������
			make_matrix(H);
			duplicate_matrix(H, H1, M, N);
			gaussian_elimination(H, H1);
			trancepose_matrix(H, Ht, M - 1, N);
			trancepose_matrix(H1, H1t, M - 1, N);
			encode(Ht, C);
			duplicate_matrix(C, C1, N - M + 1, N);//�G���[�m�F�p

			//����
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
		fprintf(fp, "%f,%lg\n", Eb_N0db, BER);//x=Eb/N0[db], y=BER���o��
	}

	//�������̊J��
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

	fclose(fp);

	return 0;
}