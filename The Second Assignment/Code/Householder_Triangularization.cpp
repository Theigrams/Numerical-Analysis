#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 20
const double eps=1e-40;
double a[N][N];
int n=10;

void input(){
	scanf("%d",&n);
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			scanf("%lf",&a[i][j]);
}

void def(){
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++){
			if(i==j)
				a[i][j]=1.52*cos(i+1.2*j);
			else if(i!=j)
				a[i][j]=sin(0.5*i+0.2*j);
		}
}

int Is_zero(double Q[N][N],int r){
	for(int i=r+2;i<=n;i++)
		if(fabs(Q[i][r])>eps)
			return 0;
	return 1;
}

int sgn(double x){
	if(x>0)
		return 1;
	return -1;
}

void Matirix_M(double A[N][N],double x[N],double y[N],double h){
	//Y=A*X/h
	for(int i=1;i<=n;i++)y[i]=0; 
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			y[i]+=A[i][j]*x[j]/h;
}

double Vector_M(double x[N],double y[N]){
	double sum=0;
	for(int i=1;i<=n;i++)
		sum+=x[i]*y[i];
	return sum;
}



void copy(double X[N][N],double Y[N][N]){
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			X[i][j]=Y[i][j];
}

void output(double A[N][N]){
	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++)
			printf("%lf	",A[i][j]);
		printf("\n");
	}
}

void Householder_Triangularization(double A[N][N]){
	double u[N],p[N],q[N],d,c,h;
	for(int r=1;r<n-1;r++){
		d=0;
		if(Is_zero(A,r))continue;
		for(int i=r+1;i<=n;i++)
			d+=A[i][r]*A[i][r];
		d=sqrt(d);
		c=-sgn(A[r+1][r])*d;
		h=c*c-c*A[r+1][r];
		for(int i=1;i<=n;i++){
			if(i<=r)
				u[i]=0;
			else if(i==(r+1))
				u[i]=A[i][r]-c;
			else if(i>r)
				u[i]=A[i][r];
		}
		memset(p,0,sizeof(p));
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				p[i]+=(A[j][i]*u[j]/h);
		Matirix_M(A,u,q,h);
		c=Vector_M(p,u)/h;
		for(int i=1;i<=n;i++)
			q[i]-=u[i]*c;
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				A[i][j]-=(q[i]*u[j]+u[i]*p[j]);
	}
}


int main(){
	double Q[N][N];
	freopen("HT.in","r",stdin);
//	freopen("HT.out","w",stdout);
	def();
	Householder_Triangularization(a);
	output(a);
	return 0;
}




