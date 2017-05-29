//在visual Studio2010的c语言环境下编译通过
#include<math.h>
#include<stdio.h>
#include<conio.h>
#include<iostream>
#define n 4
using namespace std;
static double x[n],x_[n],U[11][21],t[11][21],u[11][21],tt[8][5],uu[8][5],UU[8][5],c[9][9];
int k;
double MAX(double a[n])//求数组中的最大值
{
	int i;
	double max;
    max=fabs(a[0]);
	for(i=0;i<n;i++)
		if(fabs(a[i])>max)
			max=fabs(a[i]);		
	return(max);
}

void DooLittle(double a[n][n],double b[n])//选主元的DooLittle 分解法求线性方程组
{
	int i,j,k,t,i_k,m[n];
	double u[n][n],l[n][n],s[n],y[n];
	double u_x,l_u,max,temp;
	for(k=1;k<=n;k++)
	{ 
		for(i=k;i<=n;i++)
		{	
			l_u=0;
			for(t=1;t<=k-1;t++)  
				l_u=l_u+l[i-1][t-1]*u[t-1][k-1];
			s[i-1]=a[i-1][k-1]-l_u;
		}
		max=fabs(s[k-1]);
		i_k=k;
		m[k-1]=k;
		for(i=k;i<=n;i++)
		{
			if(fabs(s[i-1])>max)
			{
				max=fabs(s[i-1]);
				i_k=i;
			    m[k-1]=i_k;
			}
		}
		if(i_k!=k)
		{
			for(t=1;t<=k-1;t++)
			{   
				temp=l[k-1][t-1];
				l[k-1][t-1]=l[i_k-1][t-1];
				l[i_k-1][t-1]=temp;
			}
			for(t=k;t<=n;t++)
			{			
				temp=a[k-1][t-1];
				a[k-1][t-1]=a[i_k-1][t-1];
				a[i_k-1][t-1]=temp;
			}
            temp=s[k-1];
			s[k-1]=s[i_k-1];
			s[i_k-1]=temp;
		}
		u[k-1][k-1]=s[k-1];
	    for(j=k+1;j<=n;j++)
		{		
			l_u=0;
			for(t=1;t<=k-1;t++) 
				l_u=l_u+l[k-1][t-1]*u[t-1][j-1];
			u[k-1][j-1]=a[k-1][j-1]-l_u;
		}
        for(i=k+1;i<=n;i++)  
			l[i-1][k-1]=s[i-1]/u[k-1][k-1];		
	}
	for(k=1;k<=n-1;k++)	
	{
		temp=b[k-1];
		b[k-1]=b[m[k-1]-1];
		b[m[k-1]-1]=temp;
	}
	y[0]=b[0];
	for(i=2;i<=n;i++)
	{
		l_u=0;
		for(t=1;t<=i-1;t++) 
			l_u=l_u+l[i-1][t-1]*y[t-1];	
		y[i-1]=b[i-1]-l_u;
	}
	x_[n-1]=y[n-1]/u[n-1][n-1];
	for(i=n-1;i>=1;i--)
	{
		u_x=0;
		for(t=i+1;t<=n;t++)
			u_x=u_x+u[i-1][t-1]*x_[t-1];
		x_[i-1]=(y[i-1]-u_x)/u[i-1][i-1];
	}
}
void NEWTON(double x_x,double y_y)//牛顿法求解非线性方程组子程序
{ 
	double F[n],F_1[n][n];
	int i,j;
    for(i=0;i<n;i++)
	{
		x_[i]=x[i]=1;
	}
	for(i=0;i<2000;i++)
	{
		F[0]=(2.67-0.5*cos(x[0])-x[1]-x[2]-x[3]+x_x);
        F[1]=(1.07-x[0]-0.5*sin(x[1])-x[2]-x[3]+y_y);
        F[2]=(3.74-0.5*x[0]-x[1]-cos(x[2])-x[3]+x_x);
		F[3]=(0.79-x[0]-0.5*x[1]-x[2]-sin(x[3])+y_y);

        F_1[0][0]=-0.5*sin(x[0]);
 F_1[0][1]=1
; F_1[0][2]= 1;
 F_1[0][3]= 1;
F_1[1][0]=1; 
F_1[1][1]=0.5*cos(x[1]);
 F_1[1][2]= 1; 
F_1[1][3]= 1;
F_1[2][0]=0.5;
F_1[2][1]=1;
F_1[2][2]=-sin(x[2]); F_1[2][3]= 1;
F_1[3][0]=1;
F_1[3][1]=0.5;
F_1[3][2]= 1;
F_1[3][3]=cos(x[3]);
DooLittle(F_1,F);
        if(MAX(x_)/MAX(x)<=1.0e-12)
			return;		
		for(j=0;j<n;j++)
		{
			x[j]+=x_[j];
		}		
	}
	cout<<"迭代2000次没有带到精度要求"<<endl;
	return;
}
void Solve_tu()//求解对应的t、u子程序
{
	int i,j;
	double x_x[11],y_y[21];   
	for(i=0;i<11;i++)
		x_x[i]=0.08*i;	
	for(j=0;j<21;j++)
		y_y[j]=0.5+0.05*j;
    for(i=0;i<11;i++)
		for(j=0;j<21;j++)
		{      
    	    NEWTON(x_x[i],y_y[j]);
	    	t[i][j]=x[0];
			u[i][j]=x[1];
		}	
	for(i=0;i<8;i++)
		x_x[i]=0.1*(i+1);	
	for(j=0;j<5;j++)
		y_y[j]=0.5+0.2*(j+1);
    for(i=0;i<8;i++)
		for(j=0;j<5;j++)
		{    
    	    NEWTON(x_x[i],y_y[j]);
	    	tt[i][j]=x[0];
			uu[i][j]=x[1];
		}
}

void Lagrange(double *t,double *u,double *U,int la,int lb,int flag)//分片二次代数插值
{
int i,j,k,m,p,c,d,q;
double a[11][21],b[11][21],Z[6][6],temp1,temp2,L1[3],L2[3];

Z[0][0]=-0.5 ; Z[0][1]=-0.34; Z[0][2]= 0.14; Z[0][3]= 0.94; Z[0][4]= 2.06;  Z[0][5]=  3.5;
Z[1][0]=-0.42; Z[1][1]=-0.5 ; Z[1][2]=-0.26; Z[1][3]= 0.3 ; Z[1][4]= 1.18;  Z[1][5]= 2.38;
Z[2][0]=-0.18; Z[2][1]=-0.5 ; Z[2][2]=-0.5 ; Z[2][3]=-0.18; Z[2][4]= 0.46;  Z[2][5]= 1.42;
Z[3][0]= 0.22; Z[3][1]=-0.34; Z[3][2]=-0.58; Z[3][3]=-0.5 ; Z[3][4]=-0.1 ;  Z[3][5]= 0.62;
Z[4][0]= 0.78; Z[4][1]=-0.02; Z[4][2]=-0.5 ; Z[4][3]=-0.66; Z[4][4]=-0.5 ;  Z[4][5]=-0.02;
Z[5][0]= 1.5 ; Z[5][1]= 0.46; Z[5][2]=-0.26; Z[5][3]=-0.66; Z[5][4]=-0.74;  Z[5][5]=-0.5; 
for(i=0;i<la;i++)
  {
      for(j=0;j<lb;j++)
	  {       
		c=int(*(u+i*lb+j)/0.4);
		d=int((*(u+i*lb+j)-c*0.4)/(0.5*0.4));
		a[i][j]=(c+d)*0.4;		 
	
		c=int(*(t+i*lb+j)/0.2);
		d=int((*(t+i*lb+j)-c*0.2)/(0.5*0.2));
		b[i][j]=(c+d)*0.2;

		if(a[i][j]<0.4)			a[i][j]=0.4;
		if(a[i][j]>1.6)			a[i][j]=1.6;
		if(b[i][j]<0.2)			b[i][j]=0.2;
		if(b[i][j]>0.8)			b[i][j]=0.8;

		for(k=0;k<3;k++)
		{   
			temp1=1;temp2=1;
			for(m=0;m<3;m++)
				if(m!=k)
				{
					temp1*=(*(t+i*lb+j)-(b[i][j]+(m-1)*0.2))/(b[i][j]+(k-1)*0.2-(b[i][j]+(m-1)*0.2));
					temp2*=(*(u+i*lb+j)-(a[i][j]+(m-1)*0.4))/(a[i][j]+(k-1)*0.4-(a[i][j]+(m-1)*0.4));
				}
			L1[k]=temp1;
			L2[k]=temp2;
		}
		
		temp1=0;
        for(k=0;k<3;k++)
			for(m=0;m<3;m++)
			{
				p=int(b[i][j]/0.2)-1+k;
				q=int(a[i][j]/0.4)-1+m;			
				temp1+=L1[k]*L2[m]*Z[p][q];
			}
		*(U+i*lb+j)=temp1;

		if(flag)
		{
		 printf("%s%2d%s%.2f%s%2d%s%.2f%s%19.11e%s","x(",i,")=",0.08*i,", y(",j,")=",0.5+0.05*j,", ",*(U+i*lb+j),"; ");
		 if(!fmod(j+1,2))
		  	cout<<endl;
		 if(j==lb-1)
			 cout<<endl;
		}
	  }
  }
}


void Multi(double *a, double *b, double *c, int la, int lb, int lc, int r, int s, int t)
	//求r*s阶矩阵A与s*t阶矩阵B的乘积矩阵C
{
   int i, j, k;
   for (i=0; i<r; i++)
	   for (j=0; j<t; j++)
	   {
		   *(c+i*lc+j)=0;
		   for (k=0; k<s; k++)*(c+i*lc+j)+=*(a+i*la+k)*(*(b+k*lb+j));
	   }
}

double Inverse(double *a, double *b, int la, int lb, int N) //求n阶方阵A的逆矩阵B
{
   int i, j, k;
   double temp;
   for(i=0; i<N; i++)
      for(j=0; j<N; j++)
		if (i==j)
		    *(b+i*lb+j)=1;
		else
			*(b+i*lb+j)=0;
   for (k=0; k<N; k++){
      j=k;
      for (i=k+1; i<N; i++)
		  if (fabs(*(a+i*la+k))>fabs(*(a+j*la+k)))    j=i;
      if (j!=k)
	      for (i=0; i<N; i++)
		  {
	              temp=*(a+j*la+i);
	              *(a+j*la+i)=*(a+k*la+i);
	              *(a+k*la+i)=temp;
	              temp=*(b+j*lb+i);
	              *(b+j*lb+i)=*(b+k*lb+i);
	              *(b+k*lb+i)=temp;
		  }
      if (*(a+k*la+k)==0)
	      return 0;
      if ((temp=*(a+k*la+k))!=1)
	  for (i=0; i<N; i++){
		  *(a+k*la+i)/=temp;
	      *(b+k*lb+i)/=temp;
	  }
      for (i=0; i<N; i++)
		 if ((*(a+i*la+k)!=0) && (i!=k)){
	         temp=*(a+i*la+k);
	         for (j=0; j<N; j++){
				 *(a+i*la+j)-=temp*(*(a+k*la+j));
				 *(b+i*lb+j)-=temp*(*(b+k*lb+j));
			 }
		 }
   }
   return 0;
}

void solve_C()
{
     int    i,j,r,s;
     double t1[21][21], t2[21][21], t3[21][21],d[9][9],e[9][9];
     double B[11][9], B_T[9][11], G[21][9], G_T[9][21],P[11][21];
     double temp, FangCha;
     
	 for(i=0;i<9;i++)
	 {
		   for(j=0;j<11;j++)
		   {
			   B[j][i]=pow(0.08*j,i);
			   B_T[i][j]=pow(0.08*j,i);
		   }
		   for(j=0;j<21;j++)
		   {
			   G[j][i]=pow(0.5+0.05*j,i);
			   G_T[i][j]=pow(0.5+0.05*j,i);
		   }
	  }	  

	 for (k=0; k<9; k++) 
	  {
		  FangCha=0;
       	  Multi(B_T[0], B[0], t1[0], 11, 9, 21, k+1, 11, k+1);	 
		  Inverse(t1[0], c[0], 21, 9, k+1);
          Multi(e[0], c[0], d[0], 9, 9, 9, k+1, k+1, k+1);
		  Multi(c[0], B_T[0], t1[0], 9, 11, 21, k+1, k+1, 11);
		  Multi(t1[0], U[0], t2[0], 21, 21, 21, k+1, 11, 21);
		  Multi(G_T[0], G[0], t1[0], 21, 9, 21, k+1, 21, k+1);
		  Inverse(t1[0], c[0], 21, 9, k+1);
		  Multi(G[0], c[0], t3[0], 9, 9, 21, 21, k+1, k+1);
		  Multi(t2[0], t3[0], c[0], 21, 21, 9, k+1, 21, k+1);
	      for(i=0;i<11;i++)
             for(j=0;j<21;j++)
			 {
			     temp=0;
			     for(r=0;r<k+1;r++)
				   for(s=0;s<k+1;s++)
					  temp+=c[r][s]*B[i][r]*G[j][s];
			      P[i][j]=temp;      
			      FangCha+=(U[i][j]-temp)*(U[i][j]-temp);
			  }
		   printf("%s%d%s%.11e%s","k=",k,"; Sigma=",FangCa,";\n");
		  if(FangCha<=1.0e-7)
		  {
             printf("%s%d%s%.11e%s","达到精度要求时: k=",k,"; sigma=",FangCha,";\n系数c(r,s)如下：\n");
		     for(i=0;i<k+1;i++)
			 {
				 for(j=0;j<k+1;j++)
				 {
					 printf("%s%d%s%d%s%19.11e%s","c(",i,",",j,")=",c[i][j],"; ");
					 if(j==(k-1)/2)
					 cout<<endl;
				 }
				 
				 cout<<endl;

			 }
			 cout<<endl;
		  return;
		}
	}
	cout<<"经过8次拟合没有达到所需精度；"<<endl;
	return;
}

void Check()//验证子程序
{  
	int i,j,r,s;
	double B[8][9],G[5][9],temp;
  
	Lagrange( tt[0],uu[0],UU[0],8,5,0);	
    for(i=0;i<k+1;i++)
	   {
		   for(j=0;j<8;j++)
		   {
			   B[j][i]=pow(0.1*(j+1),i);
		   }
		   for(j=0;j<5;j++)
		   {
			   G[j][i]=pow(0.5+0.2*(j+1),i);
		   }
	   }	

	for(i=0;i<8;i++)
       for(j=0;j<5;j++)
		{
			 temp=0;
	         for(r=0;r<k+1;r++)
			   for(s=0;s<k+1;s++)
				  temp+=c[r][s]*B[i][r]*G[j][s];
			 printf("%s%d%s%.2f%s%d%s%.2f","x*(" , i+1, ")=" , 0.1*(i+1) , ", y*(",j+1,")=",0.5+0.2*(j+1));
             printf("%s%19.11e%s%19.11e%s" , ", f(x*,y*)=",UU[i][j],", p(x*,y*)=",temp,";\n");
	   }
}
int main()
{ 
    void NEWTON(double,double);
	void DooLittle(double **,double *);
	void Solve_tu();
	void solve_C();
	void Check();
	
	void Lagrange(double *t,double *u,double *U,int la,int lb,int );
	Solve_tu();
	Lagrange(t[0],u[0],U[0],11,21,1);
	solve_C();
Check();
getchar();
}
 

