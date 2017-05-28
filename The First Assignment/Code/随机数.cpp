#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
void srand (unsigned int n);
int n=501;
double a[600][600]; 
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
double def(int i,int j){
	int k=abs(i-j);
	if(k==0)
		return (1.64-0.024*i)*sin(0.2*i)-0.64*exp(0.1/i);
	else if(k==1)
		return 0.16;
	else if(k==2)
		return -0.064;
}

void input(){
	for(int i=1;i<=n;i++)
		for(int j=max(1,i-2);j<=min(n,i+2);j++)
			a[i][j]=def(i,j);
}
int main(){
//	freopen("input.in","r",stdin);
	freopen("input.in","w",stdout);

	for (int i=1; i<=n; i++)
		for (int j=1;j<=n;j++)
			a[i][j]=float(rand())/1000;
	memset(a,0,sizeof(a));
	input();
 	printf("%d\n", n);	

	for (int i=1; i<=n; i++)
		for (int j=1;j<=n;j++)
			printf("%lf ",a[i][j]);
	printf("\n\n\n");
	
	printf("{"); 
	for (int i=1; i<=n; i++){
		printf("{"); 
		for (int j=1;j<=n;j++){
			printf("%lf",a[i][j]);
			if(j!=n)printf(","); 
		}
		printf("}"); 
		if(i!=n)printf(","); 
	} printf("}"); 
	return 0;
}

