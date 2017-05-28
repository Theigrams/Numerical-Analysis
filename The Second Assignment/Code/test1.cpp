//#include "stdafx.h"
#include <math.h>
#include <stdio.h>
#define n 3
void output(double A[3][3]){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++)
			printf("%lf	",A[i][j]);
		printf("\n");
	}
}
int main()
{
//		freopen("QR.in","r",stdin);
	freopen("QR1.out","w",stdout);
int i,j,k,sg,r,count;
count=0;
double a[3][3],s[3],c,m,Q[3][3],R[3][3],u[3],hr,e[3],w[3],p[3],I[3][3];
a[0][0]=3;a[0][1]=1;a[0][2]=0;
a[1][0]=1;a[1][1]=2;a[1][2]=1;
a[2][0]=0;a[2][1]=1;a[2][2]=1;
while(count<10)
{	for(i=0;i<n;i++)           //Q矩阵初始化为单位阵
	{	for(j=0;j<n;j++)
		{	if(i==j)
			Q[i][j]=1;
			else
			Q[i][j]=0;
		}
	}
for(r=0;r<n-1;r++)
{	for(i=0;i<3;i++)
	{	if(i==r)
		e[i]=1;
		else
		e[i]=0;
	}
	m=0;
	for(i=0;i<r;i++)
	s[i]=0;
	for(i=r;i<n;i++)
	{	s[i]=a[i][r];
		m=m+s[i]*s[i];
	}
	if(m==0)               //如果u元素全为0，则Qr+1=Qr
	break;
	sg=(a[r][r]>0)?-1:1;      //求a[r][r]符号
	c=sg*sqrt(m);
	hr=c*c-c*a[r][r];
	printf("h=%lf\n",hr);
for(j=0;j<n;j++)
{
	for(i=0;i<n;i++)
	u[i]=s[i]-c*e[i];
}

		printf("U[]=:\n");
		for( i=0;i<n;i++)printf("%lf ",u[i]);
		printf("\n");

for(i=0;i<n;i++)
{	p[i]=0;
	for(j=0;j<n;j++)
	p[i]=p[i]+a[j][i]*u[j]/hr;
}
for(i=0;i<n;i++)
{	w[i]=0;
	for(j=0;j<n;j++)
	w[i]=w[i]+Q[i][j]*u[j]/hr;
}
	printf("--***-----\n");
		output(Q);
		printf("W[]=:\n");
		for( i=0;i<n;i++)printf("%lf ",w[i]);printf("\n");
for(i=0;i<n;i++)
{	for(j=0;j<n;j++)
	a[i][j]=a[i][j]-u[i]*p[j];
}
for(i=0;i<n;i++)
{	for(j=0;j<n;j++)
	Q[i][j]=Q[i][j]-w[i]*u[j];
}
		printf("Q[R=%d]:\n",r);
		output(Q);
}
	printf("Q[%d]:\n",count+1);
	output(Q);
	printf("R:[%d]\n",count+1);
	output(a);
for(i=0;i<n;i++)
{	for(j=0;j<n;j++)
	{	R[i][j]=a[i][j];
	}
}
for(i=0;i<n;i++)
	{	for(j=0;j<n;j++)
		R[i][j]=a[i][j];           //R=An
	}
for(i=0;i<n;i++)
{	for(j=0;j<n;j++)
	{	a[i][j]=0;
		for(k=0;k<n;k++)
		a[i][j]=a[i][j]+R[i][k]*Q[k][j];
	}
}
	count=count+1;
}
for(i=0;i<n;i++)
{	for(j=0;j<n;j++)
	printf("%lf\t",a[i][j]);
	printf("\n");
}
}

