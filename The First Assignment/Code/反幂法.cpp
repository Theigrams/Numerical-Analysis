#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 1000
double a[N][N],u[N],y[N];
const double eps=1e-10;
int n;

void input(){
	scanf("%d",&n);
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			scanf("%lf",&a[i][j]);
}

void put(double U[N]){
	for(int i=1;i<=n;i++)
		printf("%lf,",U[i]);
	printf("\n");
}

void Doolittle(){
	for(int k=1;k<=n;k++){
		for(int j=k;j<=n;j++)
			for(int t=1;t<k;t++)
				a[k][j]-=a[k][t]*a[t][j];
		for(int i=k+1;i<=n;i++){
			for(int t=1;t<k;t++)
				a[i][k]-=a[i][t]*a[t][k];
			a[i][k]/=a[k][k];
		}
	}
}

void Slove(double X[],double b[]){
	for(int i=1;i<=n;i++){
		X[i]=b[i];
		for(int t=1;t<i;t++)
			X[i]-=a[i][t]*X[t];
	}
	for(int i=n;i;i--){
		for(int t=i+1;t<=n;t++)
			X[i]-=a[i][t]*X[t];
		X[i]/=a[i][i];
	}
}

void normalization(){
	double s=0;
	for(int i=1;i<=n;i++)
		s+=u[i]*u[i];
	s=sqrt(s);
	for(int i=1;i<=n;i++)
		y[i]=u[i]/s;
}

void multiplication(double U[N]){
	memset(u,0,sizeof(u));
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			u[i]+=a[i][j]*U[j];
}

double InversePowerMethod(){
	double B=0,b=0,e=1;
	for(int i=1;i<=n;i++)
		u[i]=1;//initialization
	Doolittle();
	while(e>eps){
		
	/*	printf("U:");
		put(u);*/
		
		normalization();
		
/*		printf("y:");
		put(y);*/
		
		Slove(u,y);
		
		
		B=0;
		for(int i=1;i<=n;i++)
			B+=y[i]*u[i];
		e=fabs((B-b)/B);
		b=B;
//		printf("Lamda1=%lf\n",B);
	}
	return B;
}

int main(){
	double lamda;
	freopen("input.in","r",stdin);
	freopen("InversePowerMethod.out","w",stdout);
	memset(a,0,sizeof(a));
	memset(u,0,sizeof(u));
	input();
	lamda=InversePowerMethod();
	printf("Lamda1=%lf\n",1/lamda);
	printf("y:");
	put(y);
	return 0;
}
