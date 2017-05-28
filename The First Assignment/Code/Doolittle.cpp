#include<stdio.h>
#include<math.h>
#include<string.h>
double a[101][101];
double s[101],ss[101];
int n;
void input(){
	scanf("%d",&n);
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n+1;j++)
			scanf("%lf",&a[i][j]);
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
/*
void Slove(){
	for(int i=1;i<=n;i++)
		//y[i]=a[i][n+1]
		for(int t=1;t<i;t++)
			a[i][n+1]-=a[i][t]*a[t][n+1];
	
	for(int i=n;i;i--){
		//x[i]=a[i][n+1]
		for(int t=i+1;t<=n;t++)
			a[i][n+1]-=a[i][t]*a[t][n+1];
		a[i][n+1]/=a[i][i];
	}
}
*/
void put(){
	for(int i=1;i<=n;i++){ 
		for(int j=1;j<=n+1;j++)
		printf("%lf	",a[i][j]);
		printf("\n");
	}
}

int main(){
	freopen("input.in","r",stdin);
	freopen("Doolittle.out","w",stdout);
	memset(a,0,sizeof(a));
	input();
	for(int i=1;i<=n;i++)
		s[i]=a[i][n+1];
	Doolittle();
	Slove(ss,s);
	put();
	for(int i=1;i<=n;i++)printf("%lf,",ss[i]);
	return 0;
}
