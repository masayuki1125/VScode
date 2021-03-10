#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<complex.h> 
#define Wc 2   //列重み
#define Wr 4   //行重み
#define pi M_PI
#define LMAX 100
#define TRIAL 100
#define N 120  //列 (符号長)
#define M (N*Wc/Wr)  //行 (情報長)

//配列の初期化
void initialization_matrix(int **H){
  int i,j;

  for(i=0;i<M;i++){
    for(j=0;j<N;j++){
      H[i][j]=0;
    }
  }
}

//i行目とj行目を入れ替える
void swap_row(int **H, int i,int j){

  int k,a;
  for(k=0;k<N;k++){
    a=H[i][k];
    H[i][k]=H[j][k];
    H[j][k]=a;
  }


}

//i列目とj列目を入れ替える
void swap_column(int **H,int i,int j){


  int k,a;
  for(k=0;k<M;k++){
    a=H[k][i];
    H[k][i]=H[k][j];
    H[k][j]=a;
  }


}

//i行目にj行目を足す
void add(int **H,int i,int j){

  int k;
  
  for(k=0;k<N;k++){
    H[i][k]=H[i][k]^H[j][k];
  }

}




//ガウス消去法
void gaussian_elimination(int **H,int **H1){

  int i,j,k=0,l;
  
  //まずは上三角行列にする

  for(i=N-1;i>=(N-M);i--){
    for(j=i-(N-M);j>=0;j--){
      if(j==(i-(N-M)) && H[j][i]==0){  //対角成分が1じゃない場合、1を持つ列と入れ替える
	
	do{
	  i--;
	}while(i>=0 && H[j][i]!=1);

	if(i>=0){
	  swap_column(H,i,j+(N-M));
	  swap_column(H1,i,j+(N-M));
	  i=j+(N-M);
	}else if(i==-1 && j>1){
	  swap_row(H,j,k);
	  k++;
	  i=j+(N-M);
	  j++;
	}else if(i==-1 && j<=1) i=j+(N-M);
	
	
      
      }
      
      //対角以上が0でない場合、対角を持つ行を足しあわせる

      if((i-(N-M))!=j && H[j][i]==1){ 
	add(H,j,i-(N-M));

      }

    }

  }




/*   //行列の確認用 */
/*   for(i=0;i<M;i++){ */
/*     for(j=0;j<N;j++){ */
/*       printf("%d",H[i][j]); */
/*     } */
/*     printf("\n"); */
/*   } */
/*   printf("\n"); */



  //次に単位行列にする

  for(i=N-1;i>=(N-M);i--){
    for(j=M-1;j>=0;j--){
      if(i!=j+(N-M) && H[j][i]==1){
	add(H,j,i-(N-M));
      }
    }

  }

/*   //行列の確認用 */
/*   for(i=0;i<M;i++){ */
/*     for(j=0;j<N;j++){ */
/*       printf("%d",H[i][j]); */
/*     } */
/*     printf("\n"); */
/*   } */
/*   printf("\n"); */

}




void make_LDPC(int **H){


  int i,j,k,l;   //ループ用変数
  int a,b,c;

  //検査行列に１を挿入
  for(k=0;k<Wc;k++){
    for(i=0;i<(M/Wc);i++){
      for(j=i*Wr;j<(i+1)*Wr;j++){
	H[k*(M/Wc)+i][j]=1;
      }

    }

  }


  for(k=1;k<Wc;k++){
    for(l=0;l<5*N;l++){
      a=N*(double)rand()/RAND_MAX;
      b=N*(double)rand()/RAND_MAX;


      for(i=0;i<(M/Wc);i++){
	c=H[k*(M/Wc)+i][a];
	H[k*(M/Wc)+i][a]=H[k*(M/Wc)+i][b];
	H[k*(M/Wc)+i][b]=c;
      }
   
    }
  }

/*     //行列の確認用 */
/*     for(i=0;i<M;i++){ */
/*       for(j=0;j<N;j++){ */
/*         printf("%d",H[i][j]); */
/*       } */
/*       printf("\n"); */
/*     } */
/*     printf("\n"); */



}


void make_Gmatrix(int **H,int **G){

  int i,j,k;  //ループ用変数
  
  //行列の初期化
  for(i=0;i<(N-M+2);i++){
    for(j=0;j<N;j++){
      G[i][j]=0;
    }
  }


  //単位行列の部分
  for(i=0;i<(N-M+2);i++){
    for(j=0;j<(N-M+2);j++){
      if(i==j) G[i][j]=1;
    }
  }
  
  //行列Pの転置を代入
  for(i=0;i<(N-M+2);i++){
    for(j=N-M+2;j<N;j++){
      G[i][j]=H[j-(N-M)][i];
    }

  }
  

/*   //行列の確認用 */
/*   for(i=0;i<(N-M+2);i++){ */
/*     for(j=0;j<N;j++){ */
/*       printf("%d",G[i][j]); */
/*     } */
/*     printf("\n"); */
/*   } */
/*   printf("\n"); */

  

}



//符号化

void coding(int x[],int y[],int **G){

  int i,j,k;

  for(i=0;i<N;i++){
    y[i]=0;
  }

  for(j=0;j<N;j++){

    for(i=0;i<(N-M+2);i++){
      y[j]=y[j]^(x[i]*G[i][j]);
    }
    
    //printf("%d ",y[j]);
  }
  // printf("\n");


}


void make_signal(int x[]){

  int i;
  double u;

  //シンボルの作成
  for(i=0;i<(N-M+2);i++){
    do{
      u=(double)rand()/RAND_MAX;
    }while(u==0.0);
    if(u<0.5) x[i]=0;
    else x[i]=1;
  }


}


void BPSKmodulation(int x[],double z[]){

  int i;

  for(i=0;i<N;i++){
    if(x[i]==0) z[i]=-1;
    else if(x[i]==1) z[i]=1;

  }


}

void BPSKdemodulation(double *z,int *y){

  int i;

  for(i=0;i<N;i++){
    if(z[i]<0) y[i]=0;
    else if(z[i]>0) y[i]=1;

  }

}


void gaussiannoise(double sigma,double x[]){

  int i;
  double u1,u2;
  complex v;

  //ボックス・ミュラー
  for(i=0;i<N;i++){

  do{
    u1=(double)rand()/RAND_MAX;
  }while(u1==0.0);
  do{
    u2=(double)rand()/RAND_MAX;
  }while(u2==0.0);

  v=sqrt(-2.0*log(u1))*cos(2.0*pi*u2);


   
  x[i]=x[i]+sigma*creal(v);
  }

}




//ギャラガーのf関数
double f(double x){

  double a;
  if(x>12.5) x=12.5;
  else if(x<1E-6) x=1E-6;
  else if(x==0) x=1E-6;
  a=log((exp(x)+1)/(exp(x)-1));
  return a;

}

//sign関数
double sign(double x){

  if(x>=0) return 1;
  if(x<0)  return -1;


}




//sum-product復号法
void sum_product(int **H,int y[],int c[],double sigma){


  double **alpha;  //対数外部値比 [M][N]
  double **beta;   //対数事前値比 [M][N]
  double a,b,d,asum;
  double *temp1,*temp2,*temp3;  //予め求めておいた総和、総積
  int m,n,k,i,j;  //ループ用変数
  int l;      //反復回数
  int e;  
  int *P;  //パリティ検査に使う行列  [M]
  int **A,**B; //行列Hの1が立っている行インデックスと列インデックス  A[M][Wr] B[Wc][N]
  int **Ht;  //行列Hの転置  [N][M]

  //動的メモリ確保

  alpha=malloc(sizeof(double *) *M);
  beta=malloc(sizeof(double *) *M);

  for(m=0;m<M;m++){
    alpha[m]=malloc(sizeof(double) *N);
    beta[m]=malloc(sizeof(double) *N);
  }
    

  Ht=malloc(sizeof(int *) *N);
  for(n=0;n<N;n++){
    Ht[n]=malloc(sizeof(int) *M);
  }

  P=malloc(sizeof(int) *M);

  A=malloc(sizeof(int *) *M);
  for(m=0;m<M;m++){
    A[m]=malloc(sizeof(int) *Wr);
  }

  B=malloc(sizeof(int *) *Wc);
  for(m=0;m<Wc;m++){
    B[m]=malloc(sizeof(int) *N);
  }

  temp1=malloc(sizeof(double) *M);
  temp2=malloc(sizeof(double) *M);
  temp3=malloc(sizeof(double) *N);

  //Hを転置させる   

  for(n=0;n<N;n++){
    for(m=0;m<M;m++){
      Ht[n][m]=H[m][n];
    }
  }


  // A(m),B(n)を求める
  for(m=0;m<M;m++){
    i=0;
    for(n=0;n<N;n++){
      if(H[m][n]==1){
	A[m][i]=n;
	i++;
      }
    }
  }

  for(n=0;n<N;n++){
    i=0;
    for(m=0;m<M;m++){
      if(H[m][n]==1){
	B[i][n]=m;
	i++;
      }
    }
  }


  //step1 初期化


  for(m=0;m<M;m++){
    for(n=0;n<N;n++){
      beta[m][n]=0;
      alpha[m][n]=0;
    }
  }

  

  for(l=0;l<LMAX;l++){
    
    //A(m),B(n)を使って総和などを予め求めておく

    for(m=0;m<M;m++){
      temp1[m]=1;
      temp2[m]=0;
      for(n=0;n<Wr;n++){
	temp1[m]*=sign(2*y[A[m][n]]/(sigma*sigma)+beta[m][A[m][n]]);
	temp2[m]+=f(fabs(2*y[A[m][n]]/(sigma*sigma)+beta[m][A[m][n]]));
	
      }
    }



    //step2 行処理

    for(m=0;m<M;m++){
      for(n=0;n<Wr;n++){
	a=2*y[A[m][n]]/(sigma*sigma)+beta[m][A[m][n]];
	alpha[m][A[m][n]]=temp1[m]/sign(a)*f(temp2[m]-f(fabs(a)));
      }
	 
    }

  
    //αの総和を求めておく

    for(n=0;n<N;n++){
      temp3[n]=0;
      for(m=0;m<Wc;m++){
	temp3[n]+=alpha[B[m][n]][n];
      }
    }


    //step3 列処理

    for(n=0;n<N;n++){
      for(m=0;m<Wc;m++){
	beta[B[m][n]][n]=temp3[n]-alpha[B[m][n]][n];
      }
    }


    //step4 一時推定語の計算

    for(n=0;n<N;n++){
      if(sign(2*y[n]/(sigma*sigma)+temp3[n])==1) c[n]=1;
      else c[n]=0;

    }


    //step5 パリティ検査


    for(m=0;m<M;m++){
      P[m]=0;
    }
    e=0;


    for(m=0;m<M;m++){
      for(n=0;n<N;n++){	
	P[m]=P[m]^(c[n]*Ht[n][m]);
      }
      e=e|P[m]; 
    }

    if(e==0) break;
    

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
  for(n=0;n<N;n++){
    free(Ht[n]);
  }
  free(Ht);

  for(m=0;m<M;m++){
    free(alpha[m]);
    free(beta[m]);
  }
  free(alpha);
  free(beta);
  free(P);

  for(m=0;m<M;m++){
    free(A[m]);
  }
  free(A);

  for(m=0;m<Wc;m++){
    free(B[m]);
  }
  free(B);

  free(temp1);
  free(temp2);
  free(temp3);


}




int main(){

  int **H;  //検査行列 [M][N]
  int **H1;  //sum-product復号法で使う行列H [M][N]
  int **G;  //生成行列  [N-M+2][N]
  int *x;  //情報ビット  [N-M+2]
  int *y,*c;  //y:符号語[N] c:推定語[N]
  double *z;  //送信信号[N]
  int i,j,k,trial,sum,SUM;
  double Eb_N0,BER,ber,snr,sigma,rate=(double)Wc/(double)Wr;
  FILE *fp1;

  fp1=fopen("LDPC.txt","w");
  srand((unsigned)time(NULL));



  //動的メモリ確保

  H=malloc(sizeof(int *) *M);
  H1=malloc(sizeof(int *) *M);
  for(i=0;i<M;i++){
    H[i]=malloc(sizeof(int) *N);
    H1[i]=malloc(sizeof(int) *N);
  }

  G=malloc(sizeof(int *) *(N-M+2));
  for(i=0;i<(N-M+2);i++){
    G[i]=malloc(sizeof(int) *N);
  }

  x=malloc(sizeof(int) *(N-M+2));
  y=malloc(sizeof(int) *N);
  c=malloc(sizeof(int) *N); 
  z=malloc(sizeof(double) *N);




  initialization_matrix(H);
  make_LDPC(H);

  //検査行列Hを復号に使う行列H1にコピー
  for(i=0;i<M;i++){
    for(j=0;j<N;j++){
      H1[i][j]=H[i][j];
    }
  }


  gaussian_elimination(H,H1);
  make_Gmatrix(H,G);



  for(Eb_N0=0;Eb_N0<3;Eb_N0+=0.2){
    snr=pow(10,(Eb_N0/(double)10.0));
    sigma=sqrt((double)1.0/(double)(2*snr*rate));
    sum=0;
    for(trial=0;;trial++){
      SUM=0;

      make_signal(x);
      coding(x,y,G);
      BPSKmodulation(y,z);
      gaussiannoise(sigma,z);
      sum_product(H1,y,c,sigma);
 

      for(i=0;i<(N-M+2);i++){
	if(x[i]!=c[i]) {
	  SUM++;
	  sum++;
	}
	//printf("%d ",c[i]);
      }
      // printf("\n");

      ber=sum/(double)(N*(trial+1));
      fprintf(stderr,"Eb_N0=%f trial=%d error=%d BER=%f \r",Eb_N0,trial,SUM,ber);
      
      if(sum>300 && trial>=TRIAL) break;
      if(trial>1E+3) break;
    }


    BER=sum/(double)(trial*N);
    fprintf(fp1,"%f %f\n",Eb_N0,BER);
    printf("Eb_N0=%f trial=%d error=%d BER=%f\n",Eb_N0,trial,sum,BER);
    if(BER==0) break;
  }

  fflush(stdout);
  fflush(stderr);


  //メモリの解放

  free(x);
  free(y);
  free(z);
  free(c);

  for(i=0;i<M;i++){
    free(H[i]);
    free(H1[i]);
  }
  free(H);
  free(H1);

  for(i=0;i<N-M+2;i++){
    free(G[i]);
  }
  free(G);

}



//errorが0の時と多い時とではっきりと別れる
