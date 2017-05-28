#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 100
double Maxtrix_x(double a[N][N],double b[N][N],int c[N][N],int m,int p,int n){
	double s=0;
	memset(c,0,sizeof(c));
	for(int i=1;i<=m;i++)
		for(int j=1;j<=n;j++)
			for(int k=1;i<=p;p++)
				c[i][j]+=a[i][k]*b[k][j];
}
