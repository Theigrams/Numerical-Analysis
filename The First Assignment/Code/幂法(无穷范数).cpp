#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 600
double a[N][N],u[N],y[N];
const double eps=1e-20;
int n=501;

double mmax(double U[]){
	double m=0;
	for(int i=1;i<=n;i++)
		if(fabs(m)<fabs(U[i]))
			m=U[i];
	return m;
}

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

int normalization(){
	double s=0;
	s=mmax(u);
	for(int i=1;i<=n;i++)
		y[i]=u[i]/fabs(s);
	if(s>=0)
		return 1;
	return -1;
}

void multiplication(double U[N]){
	memset(u,0,sizeof(u));
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			u[i]+=a[i][j]*U[j];
}

double PowerMethod(){
	double B=0,b=0,e=1;
	for(int i=1;i<=n;i++)
		u[i]=15;//initialization
	while(e>eps){
		B=normalization();
		multiplication(y);
		B*=mmax(u);
		e=fabs((B-b)/B);
		b=B;
	}
	return B;
}

int main(){
	double lamda;
	freopen("input.in","r",stdin);
	freopen("11PowerMethod.out","w",stdout);
	memset(a,0,sizeof(a));
	memset(u,0,sizeof(u));
	input();
	lamda=PowerMethod();
	printf("Lamda1=%lf\n",lamda-10.7);
	printf("y:");
	put(y);
	return 0;
}
