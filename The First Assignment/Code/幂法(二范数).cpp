#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 1000
double a[N][N],u[N],y[N];
const double eps=0.0000001;
int n;
double abs(double x){
	if(x>0)
		return x;
	return -x;
}

void input(){
	scanf("%d",&n);
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			scanf("%lf",&a[i][j]);
}

void put(double U[N]){
	for(int i=1;i<=n;i++)
		printf("%5e,",U[i]);
	printf("\n");
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

double PowerMethod(){
	double B=0,b=0,e=1;
	u[1]=1;//initialization
	while(e>eps){
		normalization();
		multiplication(y);
		B=0;
		for(int i=1;i<=n;i++)
			B+=y[i]*u[i];
		e=abs((B-b)/B);
		b=B;
	}
	return B;
}

int main(){
	double lamda;
	freopen("input.in","r",stdin);
	freopen("PowerMethod.out","w",stdout);
	memset(a,0,sizeof(a));
	memset(u,0,sizeof(u));
	input();
	lamda=PowerMethod();
	printf("Lamda1=%lf\n",lamda);
	printf("y:");
	put(y);
	return 0;
}
