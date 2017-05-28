#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 10
double a[N][N];
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

void output(double A[N][N]){printf("{");
	for(int i=1;i<=n;i++){printf("{");
		for(int j=1;j<=n;j++)
			{printf("%lf ",A[i][j]);if(j!=n)printf(",");}
			printf("}");
			if(i!=n)printf(",");
		printf("\n");
	}printf("}");
}
int main(){
	double Q[N][N];
	freopen("Works.in","r",stdin);
	freopen("oy.out","w",stdout);
	def();

	output(a);
	return 0;
}

