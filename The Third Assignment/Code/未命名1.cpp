#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 100
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
	}
	printf("wrong!\n");
}

double Interpolation(double t,double u){
	double h=0.2,tau=0.4;
	double z[N][N],sum=0,p;
	int i,j;
z[0][0]=-0.5 ; z[0][1]=-0.34; z[0][2]= 0.14; z[0][3]= 0.94; z[0][4]= 2.06;  z[0][5]= 3.5;
z[1][0]=-0.42; z[1][1]=-0.5 ; z[1][2]=-0.26; z[1][3]= 0.3 ; z[1][4]= 1.18;  z[1][5]= 2.38;
z[2][0]=-0.18; z[2][1]=-0.5 ; z[2][2]=-0.5 ; z[2][3]=-0.18; z[2][4]= 0.46;  z[2][5]= 1.42;
z[3][0]= 0.22; z[3][1]=-0.34; z[3][2]=-0.58; z[3][3]=-0.5 ; z[3][4]=-0.1 ;  z[3][5]= 0.62;
z[4][0]= 0.78; z[4][1]=-0.02; z[4][2]=-0.5 ; z[4][3]=-0.66; z[4][4]=-0.5 ;  z[4][5]=-0.02;
z[5][0]= 1.5 ; z[5][1]= 0.46; z[5][2]=-0.26; z[5][3]=-0.66; z[5][4]=-0.74;  z[5][5]=-0.5; 
	if(t<0.3)i=1;
	else if(t>0.7)i=4;
	else if(0.3<=t && t<=0.7)i=int(t/h+0.5);
	if(u<0.6)j=1;
	else if(u>1.4)j=4;
	else if(0.6<=u && u<=1.4)j=int(u/tau+0.5);

	for(int k=i-1;k<=i+1;k++)
		for(int r=j-1;r<=j+1;r++){
			p=1;
			for(int m=i-1;m<=i+1;m++){
				if(m==k)continue;
				p*=(t-h*m)/(h*k-h*m);
			}
			for(int n=j-1;n<=j+1;n++){
				if(n==r)continue;
				p*=(u-tau*n)/(tau*r-tau*n);
			}
			p*=z[k][r];
			sum+=p;
		}
	return sum;
}

double Surface_Fitting(double x[N],double y[N],double z[N][N],double C[N],[N],int k){
	double B[N][N],G[N][N],I[N][N],J[N][N],sigma;
	for(int i=1;i<=11;i++)
		for(int j=1;j<=k+1;j++)
			B[i][j]=pow(x[i-1],j-1);
	for(int i=1;i<=21;i++)
		for(int j=1;j<=k+1;j++)
			G[i][j]=pow(y[i-1],j-1);
	Transpose(B,I,11,k+1);
	memset(J,0,sizeof(J));
	Maxtrix_x(I,B,J,k+1,11,k+1);
	memset(C,0,sizeof(C));
	Inverse(J,C,k+1);
	memset(J,0,sizeof(J));
	Maxtrix_x(C,I,J,k+1,k+1,11);
	memset(I,0,sizeof(I));
	Maxtrix_x(J,z,I,k+1,11,21);
	memset(J,0,sizeof(J));
	Maxtrix_x(I,G,J,k+1,21,k+1);
	memset(I,0,sizeof(I));
	Transpose(G,I,21,k+1);
	memset(C,0,sizeof(C));
	Maxtrix_x(I,G,C,k+1,21,k+1);
	memset(I,0,sizeof(I));
	Inverse(C,I,k+1);
	memset(C,0,sizeof(C));
	Maxtrix_x(J,I,C,k+1,k+1,k+1);
	
	
}
void output(){
	for(int i=1;i<=4;i++)
	printf("X[%d]=%lf\n",i,X[i]);
}

int	main(){
//	freopen("input.in","r",stdin);
//	freopen("Newton.out","w",stdout);
//	memset(a,0,sizeof(a));
	double x[N],y[N],z[N][N],C[N][N];
	for(int i=0;i<=10;i++)
		for(int j=0;j<=20;j++){
			x[i]=0.08*i;
			y[j]=0.5+0.05*j;
			for(int k=1;k<=4;k++)X[k]=1;
			Newton(x[i],y[j]);
			z[i][j]=Interpolation(X[1],X[2]);
		//	printf("f(X%d,Y%d)=%lf\n",i,j,z);
		}
	for(int k=0;k<=10;k++){
		Surface_Fitting(x,y,z,C);
		
	}
	return 0;
}

