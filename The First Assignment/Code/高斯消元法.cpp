#include<stdio.h>
#include<math.h>
#include<string.h>
double a[101][101],b[101][101],s[101];
int n,flag=0;
void input(){
	scanf("%d",&n);
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n+1;j++)
			scanf("%lf",&a[i][j]);
}

void put(){
	for(int i=1;i<=n;i++){ 
		for(int j=1;j<=n+1;j++)
		printf("%lf	",b[i][j]);
		printf("\n\n\n");}
}


void guass(int x){
	double m;
		if(x==n)return ;
		for(int i=x+1;i<=n;i++){
				if(b[x][x]==0)continue;
				m=b[i][x]/b[x][x];
					for(int j=x;j<=n+1;j++)
						b[i][j]-=m*b[x][j];
		}
	guass(x+1);
}


void slove(){
	double sum;
	for(int i=n;i>=1;i--){
		sum=0;	
		for(int j=i+1;j<=n;j++) 
			sum+=b[i][j]*s[j];
		if(b[i][i])
			s[i]=(b[i][n+1]-sum)/b[i][i];
		else
			s[i]=0;
	}
}


void check(){
	double sum;
	for(int i=1;i<=n;i++){
		sum=0;
		for(int j=1;j<=n;j++)
			sum+=s[j]*a[i][j];
		if(abs(sum -a[i][n+1])>0.01){
			printf("n=%d ,sum=%lf but ans=%lf\n",i,sum,a[i][n+1]);
			return;
		}
	}
	flag=1;
}

void output(){
	for(int i=1;i<=n;i++)
	printf("X[%d]=%lf\n",i,s[i]);
}

int	main(){
	freopen("input.in","r",stdin);
	freopen("guass1.out","w",stdout);
	memset(a,0,sizeof(a));
	memset(s,0,sizeof(s));
	input();
	memcpy(b,a,sizeof(a));
	guass(1);
	slove();
	check();
	if(flag==1)
		output();
	else
		printf("There is no answer!!!");
	return 0;
}

