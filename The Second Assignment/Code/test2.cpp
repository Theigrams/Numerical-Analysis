#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 20
int n=10;
void def(double a[N][N],double b[N]){
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++){
			if(i==j)
				a[i][j]=1.52*cos(i+1.2*j);
			else if(i!=j)
				a[i][j]=sin(0.5*i+0.2*j);
		}
	for(int i=1;i<=n;i++)b[i]=1;
}

void clear(double a[N][N],double b[N]){
	memset(a,0,sizeof(a));
	memset(b,0,sizeof(b));
}
void Q(){
	double a[N][N],b[N];
	def(a,b);
	clear(a,b);
/*	memset(a,0,sizeof(a));
	memset(b,0,sizeof(b));*/
	for(int i=1;i<=n;i++){
		for(int j=1;j<=n;j++)
			printf("%lf	",a[i][j]);
		printf("\n");
	}
	for(int i=1;i<=n;i++)printf("%lf	",b[i]);
}

int main(){
	Q();
	return 0;
}
