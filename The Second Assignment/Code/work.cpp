#include<stdio.h>
#include<math.h>
#include<string.h>
#define N 20
const double eps=1e-12;
double a[N][N],B[N][N],C[N][N];
int n=10;

typedef struct{
//定义复数结构体
    double Re;
    double Im;
}ComplexNumber;

void zeroMat(double a[N][N]){
    for(int i= 1; i<=10; i++) {
        for(int j = 1; j<=10; j++) {
            if (fabs(a[i][j])<eps)
                a[i][j] = 0;
        }
    }
}

void def(){
	//初始化A数组
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++){
			if(i==j)
				a[i][j]=1.52*cos(i+1.2*j);
			else if(i!=j)
				a[i][j]=sin(0.5*i+0.2*j);
		}
}

void eye(double Q[N][N]){
	//定义单位阵 
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			Q[i][j]=(i==j);
}

int Is_zero(double Q[N][N],int r,int k){
	//K>0表示在QR中调用，K<0表示在 DQR中调用 
	int kk=n;//更改数组大小 
	if(k<0){kk=-k;k=1;}
	for(int i=r+k;i<=kk;i++)
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

void MatirixM_M(double X[N][N],double Y[N][N]){
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			for(int k=1;k<=n;k++)
				B[i][j]+=X[i][k]*Y[k][j];
}

void copy(double X[N][N],double Y[N][N]){
	for(int i=1;i<=n;i++)
		for(int j=1;j<=n;j++)
			X[i][j]=Y[i][j];
}

void output(double A[N][N]){
	for(int i=1;i<=n;i++){
		for(int j=6;j<=10;j++){ 
			printf("%.12e",A[i][j]);
			if(j!=10)printf(" & ");
			else if(j==10)printf("\\\\");
		}
		//printf("\n\n");
	}
}

void QR(double A[N][N],double Q[N][N]){
	double u[N],w[N],p[N],d,c,h;
	eye(Q);
	for(int r=1;r<n;r++){
		d=0;
		if(Is_zero(A,r,1))continue;
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

void Householder_Triangularization(double A[N][N]){
	double u[N],p[N],q[N],d,c,h;
	for(int r=1;r<n-1;r++){
		d=0;
		if(Is_zero(A,r,2))continue;
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


void select(int k,double b[N][N]){
	double max,c;
	int m=k;
	max=b[k][k];
	for(int i=k+1;i<=n;i++)
		if(fabs(b[i][k])>fabs(max)){
			max=b[i][k];
			m=i;
		}
	if(m!=k)
		for(int i=k;i<=n+1;i++){
		c=b[m][i];
		b[m][i]=b[k][i];
		b[k][i]=c;
		}
}

//Gauss消元法
void gauss(double lambda)
{
    double X[N];
    double m,sum,t;
    def();
    memset(X,0,sizeof(X));
    for(int i=1;i<=n;i++){
        a[i][i]-=lambda;
        a[i][n+1]=0;
    }
	for(int k=1;k<n;k++){
		select(k,a);
		for(int i=k+1;i<=n;i++){
			m=a[i][k]/a[k][k];
		for(int j=k;j<=n+1;j++)
			a[i][j]-=m*a[k][j];
		}	
	}
	X[n]=1;
	for(int i=n-1;i>=1;i--){
		sum=0;	
		for(int j=i+1;j<=n;j++) 
			sum+=a[i][j]*X[j];
		if(fabs(a[i][i])>eps)
			X[i]=(a[i][n+1]-sum)/a[i][i];
		else
			X[i]=0;
	}
	t=0;
	for(int i=1;i<=n;i++)
		t+=X[i]*X[i];
	t=sqrt(t);
	printf("Eigenvector=(");
	for(int i=1;i<=n;i++){
        printf("%lf", X[i]/t);
        if(i!=n)printf(",\\quad ");}
    printf(")\n\n");
}

void DQR(double A[N][N]){
	double Q[N][N],M[N][N],s,t,re,im,a,b,c,d;
	double u[N],v[N],p[N],q[N],h,det;
	ComplexNumber L[N]; 
	int m=n,LL=1000,r=1;
	for(int k=1;k<=100;k++){
		if (m == 1) {
            L[r].Re = A[m][m];
            L[r].Im = 0; break;
        }
        else if(m<=0)break; 
		if(fabs(A[m][m-1])<=eps){
			L[r].Re = A[m][m];
            L[r].Im = 0;
            m--;r++;continue;
		}
		a=A[m-1][m-1];
		b=A[m-1][m];
		c=A[m][m-1];
		d=A[m][m];
		re=(a+d)/2;
		det=(a-d)*(a-d)+4*b*c;
		
		if((m==2)||(fabs(A[m-1][m-2])<=eps)){
			if(det>0){
				L[r].Re = re+sqrt(det)/2;
            	L[r].Im = 0;
            	L[r+1].Re = re-sqrt(det)/2;
            	L[r+1].Im = 0;
            }
            else if(det<0){
            	L[r].Re = re;
            	L[r].Im=sqrt(fabs(det))/2;
            	L[r+1].Re = re;
            	L[r+1].Im = -L[r].Im;
			}
            m-=2;
			r+=2;
			continue;
		}
		if(k==LL)break;
		s=a+d;
		t=a*d-c*b;
		memset(M,0,sizeof(M));
		for(int i=1;i<=m;i++)
			for(int j=1;j<=m;j++){
				for(int l=1;l<=m;l++)
					M[i][j]+=A[i][l]*A[l][j];
				M[i][j] -= s*A[i][j];
				M[i][j] += t*(i==j);
			}

		for(int r=1;r<m;r++){
			d=0;
			if(Is_zero(M,r,-m))continue;
			for(int i=r;i<=m;i++)
				d+=M[i][r]*M[i][r];
			d=sqrt(d);
			c=-sgn(M[r][r])*d;
			h=c*c-c*M[r][r];
			for(int i=1;i<=m;i++){
				if(i<r)
					u[i]=0;
				else if(i==r)
					u[i]=M[r][r]-c;
				else if(i>r)
					u[i]=M[i][r];
			}
			memset(p,0,sizeof(p));
			memset(q,0,sizeof(q));
			memset(v,0,sizeof(v));
			for(int i=1;i<=m;i++)
				for(int j=1;j<=m;j++)
					v[i]+=(M[j][i]*u[j]/h);	
			for(int i=1;i<=m;i++)
				for(int j=1;j<=m;j++){
					M[i][j]-=(u[i]*v[j]);
					p[i]+=(A[j][i]*u[j]/h);
					q[i]+=(A[i][j]*u[j]/h);
				}

			c=0;
			for(int i=1;i<=m;i++)
				c+=p[i]*u[i]/h;
			for(int i=1;i<=m;i++)
				q[i]-=u[i]*c;
			for(int i=1;i<=m;i++)
				for(int j=1;j<=m;j++)
					A[i][j]-=(q[i]*u[j]+u[i]*p[j]);
			zeroMat(A);
		}

	}
	zeroMat(A);
 	for(int r = 1; r<=10; r++){
        printf("\n\n");
        if (L[r].Im == 0) {
            printf("$\\lambda_{%d}= %.12e$ \n\n", r , L[r].Re); 
            gauss(L[r].Re);
        }
        else if(L[r].Im>0){
            printf("$\\lambda_{%d}= %.12e$+ i*%.12e\n\n", r , L[r].Re, L[r].Im);
        }
        else if(L[r].Im<0){
            printf("$\\lambda_{%d}= %.12e$- i*%.12e\n\n", r , L[r].Re, -L[r].Im);
        }
    }	
}

int main(){
	double Q[N][N];
	freopen("Works.in","r",stdin);
	freopen("Works.out","w",stdout);
	def();
	Householder_Triangularization(a);
	printf("A_n-1:\n");
	output(a);
	QR(a,Q);
	printf("\nR:\n\n");
	output(a);
	printf("\nQ:\n\n");
	output(Q);
	memset(B,0,sizeof(B));
	MatirixM_M(a,Q);
	printf("\nR*Q:\n");
	output(B);
	def();
	DQR(a);
	return 0;
}




