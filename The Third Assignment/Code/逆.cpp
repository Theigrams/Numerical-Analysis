#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 100
double a[N][N]; 
int n,flag=0;
void input(){
	scanf("%d",&n);
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			scanf("%lf",&a[i][j]);
}

void put(double b[N][N]){
	for(int i=1;i<=n;i++){ 
		for(int j=1;j<=n;j++)
		printf("%lf	",b[i][j]);
		printf("\n");
	}
}


void Inverse(double A[N][N],double B[N][N]){
	double m;
	memset(B,0,sizeof(B));
	for(int k=1;k<=n;k++)B[k][k]=1;
	for(int k=1;k<n;k++){
		for(int i=k+1;i<=n;i++){
			m=A[i][k]/A[k][k];
			for(int j=1;j<=n;j++){
				A[i][j]-=m*A[k][j];
				B[i][j]-=m*B[k][j];
			}
		}
	} 
	for(int k=n;k;k--){
		for(int i=n;i>k;i--){
			m=a[k][i];
			for(int j=1;j<=n;j++){
				A[k][j]-=m*A[i][j];
				B[k][j]-=m*B[i][j];
			}
		}
		m=A[k][k];
		for(int j=1;j<=n;j++){
			A[k][j]/=m;
			B[k][j]/=m;
		}
	}
}



int	main(){
	double b[N][N];
	freopen("input.in","r",stdin);
	freopen("Inverse.out","w",stdout);
 
	memset(a,0,sizeof(a));
	input();
	Inverse(a,b);
	put(a);
	put(b);
	return 0;
}

