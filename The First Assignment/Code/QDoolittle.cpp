#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 200
double a[101][101];
int M[101];
int n;
void input(){
	scanf("%d",&n);
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n+1;j++)
			scanf("%lf",&a[i][j]);
	for(int i=1;i<=n;i++)
		M[i]=i;
}
void select(int k,double A[N][N]){
	double max,c,sum;
	int m=k;
	max=0;
	for(int i=k;i<=n;i++){
		sum=A[i][k];
		for(int t=1;t<=k-1;t++)
			sum-=A[i][t]*A[t][k];
		if(abs(sum)>abs(max)){
			max=sum;
			m=i;
		}
	}
	if(m!=k)
		for(int t=1;t<=n;t++){
			c=A[m][t];
			A[m][t]=A[k][t];
			A[k][t]=c;
		}
	A[k][k]=max;
	M[k]=m;
}

void Qb(double A[N][N]){
	double t;
	for(int k=1;k<=n;k++){
		t=A[M[k]][n+1];
		A[M[k]][n+1]=A[k][n+1];
		A[k][n+1]=t;	
	}
}

void Doolittle(double A[N][N]){
	for(int k=1;k<=n;k++){
		select(k,A);
		for(int j=k+1;j<=n;j++)
			for(int t=1;t<k;t++)
				A[k][j]-=A[k][t]*A[t][j];
		for(int i=k+1;i<=n;i++){
			for(int t=1;t<k;t++)
				A[i][k]-=A[i][t]*a[t][k];
			A[i][k]/=A[k][k];
		}
	}
}

void Slove(double X[],double b[]){
	for(int i=1;i<=n;i++){
		X[i]=b[i];
		for(int t=1;t<i;t++)
			X[i]-=A[i][t]*X[t];
	}
	for(int i=n;i;i--){
		for(int t=i+1;t<=n;t++)
			X[i]-=A[i][t]*X[t];
		X[i]/=A[i][i];
	}
}

void put(){
	for(int i=1;i<=n;i++){ 
		for(int j=1;j<=n+1;j++)
		printf("%lf	",a[i][j]);
		printf("\n");
	}
}

int main(){
	freopen("input.in","r",stdin);
	freopen("QDoolittle.out","w",stdout);
	memset(a,0,sizeof(a));
	input();
	Doolittle();
	Qb();
	Slove();
	put();
	return 0;
}
