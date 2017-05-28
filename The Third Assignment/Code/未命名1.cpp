#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 10
double X[N];//t,u,v,w 
void F(double x,double y,double Fx[N]){
	Fx[1]=0.5*cos(X[1])+X[2]+X[3]+X[4]-x-2.67;
	Fx[2]=X[1]+0.5*sin(X[2])+X[3]+X[4]-y-1.07;
	Fx[3]=0.5*(X[1])+X[2]+cos(X[3])+X[4]-x-3.74;
	Fx[4]=X[1]+0.5*X[2]+X[3]+sin(X[4])-y-0.79;
}

void JF(double x,double y,double A[N][N]){
	A[1][1]=-0.5*sin(X[1]);	A[1][2]=1;	A[1][3]=1;	A[1][4]=1;
	A[2][1]=1;	A[2][2]=0.5*cos(X[2]);	A[2][3]=1;	A[2][4]=1;
	A[3][1]=0.5;A[3][2]=1;	A[3][3]=-sin(X[3]);		A[3][4]=1;
	A[4][1]=1;	A[4][2]=0.5;	A[4][3]=1;	A[4][4]=cos(X[4]);
}

double Maxtrix_x(double a[N][N],double b[N][N],int c[N][N],int m,int p,int n){
	double s=0;
	memset(c,0,sizeof(c));
	for(int i=1;i<=m;i++)
		for(int j=1;j<=n;j++)
			for(int k=1;i<=p;p++)
				c[i][j]+=a[i][k]*b[k][j];
}

void Inverse(double C[N][N],double B[N][N],double n){
	double m,A[N][N];
	memset(B,0,sizeof(B));
	memset(A,0,sizeof(A));
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			A[i][j]=C[i][j];
	for(int k=1;k<=n;k++)B[k][k]=1;
	for(int k=1;k<n;k++){
		for(int i=k+1;i<=n;i++){
			m=A[i][k]/A[k][k];
			for(int j=1;j<=n;j++){
				A[i][j]-=m*A[k][j];
				B[i][j]-=m*B[k][j];
			}
		}
	} 
	for(int k=n;k;k--){
		for(int i=n;i>k;i--){
			m=A[k][i];
			for(int j=1;j<=n;j++){
				A[k][j]-=m*A[i][j];
				B[k][j]-=m*B[i][j];
			}
		}
		m=A[k][k];
		for(int j=1;j<=n;j++){
			A[k][j]/=m;
			B[k][j]/=m;
		}
	}
}

void Transpose(double A[N][N],double B[N][N],int n,int m){
	memset(B,0,sizeof(B));
	for(int i=1;i<=n;i++)
		for(int j=1;j<=m;j++)
			B[i][j]=A[j][i];
}

double MaxX(double x[N],int n){
	double max=0;
	for(int i=1;i<=n;i++)
		if(fabs(x[i])>max)max=fabs(x[i]);
}

void Newton(double x,double y){
	int n=4,k;
	double eps=1e-12,max,A[N][N],Fx[N],B[N][N],dX[N];
	for(k=1;k<=10000;k++){
		memset(Fx,0,sizeof(Fx));
		memset(B,0,sizeof(B));
		memset(dX,0,sizeof(dX));
		F(x, y, Fx);
		JF(x, y, A);
		Inverse(A,B,n);
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				dX[i]+=Fx[j]*B[i][j];
		max=0;
		for(int i=1;i<=n;i++)
			if(fabs(X[i])>max)max=fabs(X[i]);
		if((MaxX(dX,n)/max)<eps)
			return ;
		for(int i=1;i<=n;i++)X[i]-=dX[i];

for(int i=1;i<=4;i++)
printf("dX[%d]=%lf,",i,dX[i]);
printf("\n");
for(int i=1;i<=4;i++)
printf("X[%d]=%lf,",i,X[i]);
printf("\n");
	}
	if(k>=99)printf("wrong!\n");
}
 
void output(){
	for(int i=1;i<=4;i++)
	printf("X[%d]=%lf\n",i,X[i]);
}

int	main(){
//	freopen("input.in","r",stdin);
//	freopen("Newton.out","w",stdout);
//	memset(a,0,sizeof(a));
	double x=0.08,y=0.55;
	for(int i=1;i<=4;i++)X[i]=1;
	Newton(x,y);
	output();
	return 0;
}

