#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 510
double u[N],y[N],c[10][N];
const double eps=1e-12;
int n=501;

double A(int i,int j){
	int k=i-j+3;
	if(k<=0)
		return 0;
	return c[k][j];
}

double mmax(double U[]){
	double m=0;
	for(int i=1;i<=n;i++)
		if(fabs(m)<fabs(U[i]))
			m=U[i];
	return m;
}

int max(int x,int y){
	if(x>y)
		return x;
	return y;
}

int min(int x,int y){
	if(x<y)
		return x;
	return y;
}

void def(double s){
	for(int j=1;j<=n;j++){
		c[3][j]=(1.64-0.024*j)*sin(0.2*j)-0.64*exp(0.1/j)-s;
		c[2][j]=c[4][j]=0.16;
		c[1][j]=c[5][j]=-0.064;
	}
	
}


void put(double U[N]){
	for(int i=1;i<=n;i++)
		printf("%lf,",U[i]);
	printf("\n");
}

void Doolittle(){
	for(int k=1;k<=n;k++){
		for(int j=k;j<=min(k+2,n);j++)
			for(int t=max(1,j-2);t<k;t++)
				c[k-j+3][j]-=A(k,t)*A(t,j);
		for(int i=k+1;i<=min(k+2,n);i++){
			for(int t=max(1,i-2);t<k;t++)
				c[i-k+3][k]-=A(i,t)*A(t,k);
			c[i-k+3][k]/=A(k,k);
		}
	}
}

void Slove(double X[],double b[]){
	for(int i=1;i<=n;i++){
		X[i]=b[i];
		for(int t=max(1,i-2);t<i;t++)
			X[i]-=A(i,t)*X[t];
	}
	for(int i=n;i;i--){
		for(int t=i+1;t<=min(i+2,n);t++)
			X[i]-=A(i,t)*X[t];
		X[i]/=A(i,i);
	}
}


void Normalization(){
	double s=0;
	for(int i=1;i<=n;i++)
		s+=u[i]*u[i];
	s=sqrt(s);
	for(int i=1;i<=n;i++)
		y[i]=u[i]/s;
}

void Multiplication(double U[N]){
	memset(u,0,sizeof(u));
	for(int i=1;i<=n;i++)
		for(int j=max(1,i-2);j<=min(n,i+2);j++)
			u[i]+=A(i,j)*U[j];
}

double PowerMethod(){
	double B=0,b=0,e=1;
	for(int i=1;i<=n;i++)
		u[i]=1;
	while(e>eps){
		Normalization();
		Multiplication(y);
		B=0;
		for(int i=1;i<=n;i++)
			B+=y[i]*u[i];
		e=fabs((B-b)/B);
		b=B;
	}
	return B;
}

double InversePowerMethod(){
	double B=0,b=0,e=1;
	for(int i=1;i<=n;i++)
		u[i]=1;
	Doolittle();
	while(e>eps){		
		Normalization();	
		Slove(u,y);		
		B=0;
		for(int i=1;i<=n;i++)
			B+=y[i]*u[i];
		e=fabs((B-b)/B);
		b=B;
	}
	return B;
}

int main(){
	double lambda_1,lambda_S,lambda_501,uk,det=1;
//	freopen("input.in","r",stdin);
//	freopen("Works.out","w",stdout);
	memset(c,0,sizeof(c));
	memset(u,0,sizeof(u));
	def(0);
	lambda_1=PowerMethod();
	lambda_S=1/InversePowerMethod();
	for(int i=1;i<=n;i++)
		det*=A(i,i);
	def(lambda_1);
	lambda_501=PowerMethod()+lambda_1;
	printf("(1):\n");
	printf("Lambda_1=%.12e\n",lambda_1);
	printf("Lambda_S=%.12e\n",lambda_S);
	printf("Lambda_501=%.12e\n",lambda_501);
	printf("-------------------------------\n(2):\n");
	for(int k=1;k<=39;k++){
		uk=lambda_1+k*(lambda_501-lambda_1)/40;
		def(uk);
		printf("Lambda_i%d=%.12e\n",k,1/InversePowerMethod()+uk);
	}
	printf("-------------------------------\n(3):\n");
	printf("Cond(A)=%.12e\n",fabs(lambda_1/lambda_S));
	printf("Det(A)=%.12e\n",fabs(det));
	return 0;
}
