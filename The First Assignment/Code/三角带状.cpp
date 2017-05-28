#include<stdio.h>
#include<math.h>
#include<string.h>
double a[101][101];
int n,s,r;
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

void input(){
	scanf("%d%d%d",&n,&s,&r);
	for(int i=1;i<=n;i++)
		for(int j=max(1,i-r);j<=min(n,1+s);j++)
			scanf("%lf",&a[i][j]);
}

void Doolittle(){
	for(int k=1;k<=n;k++){
		for(int j=k;j<=min(k+s,n);j++) 
			for(int t=max(max(1,k-r),j-s);t<k;t++)
				a[k][j]-=a[k][t]*a[t][j];
		for(int i=k+1;i<=min(k+r,n);i++){
			for(int t=max(max(1,i-r),k-s);t<k;t++)
				a[i][k]-=a[i][t]*a[t][k];
			a[i][k]/=a[k][k];
		}
	}
}

void Slove(){
	for(int i=1;i<=n;i++)
		for(int t=max(1,i-r);t<i;t++)
			a[i][n+1]-=a[i][t]*a[t][n+1];
	
	for(int i=n;i;i--){
		for(int t=i+1;t<=min(i+s,n);t++)
			a[i][n+1]-=a[i][t]*a[t][n+1];
		a[i][n+1]/=a[i][i];
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
	freopen("banded_linear_equations.out","w",stdout);
	memset(a,0,sizeof(a));
	input();
	Doolittle();
	Slove();
	put();
	return 0;
}
