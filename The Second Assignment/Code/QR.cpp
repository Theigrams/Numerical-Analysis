#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 20
const double eps=1e-20;
double a[N][N],b[N][N];
int n=10;

void def(){
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++){
			if(i==j)
				a[i][j]=1.52*cos(i+1.2*j);
			else if(i!=j)
				a[i][j]=sin(0.5*i+0.2*j);
		}
}

void input(){
	scanf("%d",&n);
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			scanf("%lf",&a[i][j]);
}

void eye(double Q[N][N]){
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			Q[i][j]=0;
	for(int i=1;i<=n;i++)
		Q[i][i]=1;
}

int Is_zero(double Q[N][N],int r){
	for(int i=r+1;i<=n;i++)
		if(fabs(Q[i][r])>eps)
			return 0;
	return 1;
}

int sgn(double x){
	if(x>0)
		return 1;
	return -1;
}
 
double Vector_M(double x[N],double y[N]){
	double sum=0;
	for(int i=1;i<=n;i++)
		sum+=x[i]*y[i];
	return sum;
}

void Matirix_M(double A[N][N],double x[N],double y[N],double h){
	for(int i=1;i<=n;i++)y[i]=0; 
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			y[i]+=A[i][j]*x[j]/h;
}

void MatirixM_M(double X[N][N],double Y[N][N]){

	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			for(int k=1;k<=n;k++)
				b[i][j]+=X[i][k]*Y[k][j];
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
 
void QR(double A[N][N],double Q[N][N]){
	double u[N],w[N],p[N],d,c,h;
	eye(Q);
	for(int r=1;r<n;r++){
		d=0;
		if(Is_zero(A,r))continue;
		for(int i=r;i<=n;i++)
			d+=A[i][r]*A[i][r];
		d=sqrt(d);
		c=-sgn(A[r][r])*d;
		h=c*c-c*A[r][r];
		for(int i=1;i<=n;i++){
			if(i<r)
				u[i]=0;
			else if(i==r)
				u[i]=A[r][r]-c;
			else if(i>r)
				u[i]=A[i][r];
		}
		memset(w,0,sizeof(w));
		Matirix_M(Q,u,w,1);
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				Q[i][j]-=(w[i]*u[j]/h);
		memset(p,0,sizeof(p));
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				p[i]+=(A[j][i]*u[j]/h);
		for(int i=1;i<=n;i++)
			for(int j=1;j<=n;j++)
				A[i][j]-=(u[i]*p[j]);
	}
}



int main(){
	double Q[N][N];
//	freopen("QR.in","r",stdin);
//	freopen("QR.out","w",stdout);
	def();
	for(int i=1;i<=1;i++){
		QR(a,Q);
	//	output(a);
		memset(b,0,sizeof(b));
		MatirixM_M(a,Q);
	printf("R:\n");
	output(a);
		copy(a,b);
	}

	printf("Q:\n");
	output(Q);
		printf("RQ:\n");
	output(a);
	return 0;
}




