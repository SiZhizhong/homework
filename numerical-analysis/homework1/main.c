#include <stdio.h>
#include <stdlib.h>
#include<math.h>

#define SIZE 501

//注明：由于题目中给出的sin(0.2i)中，没有注明是弧度还是角度，本算法全部按照弧度计算。


/***************************************************************
幂法部分
*************************************************************/

// 向量除以相邻的模长，对应书中的第一个和第二个公式,将向量归一化
//v 被归一化的向量   nomalized_v  归一化之后的结果
void vector_normalize(double v[SIZE],double nomalized_v[SIZE])
{
    double length=0;
    int i=0;

    for(i=0;i<SIZE;i++)
    {
        length=length+v[i]*v[i];
    }


    for(i=0;i<SIZE;i++)
        nomalized_v[i]=v[i]/sqrt(length);
}

//A矩阵乘以向量v,由于要求A的非0数据都不保存，故单独定义这种乘法，其中非0元素均不保存
//v  A*v中的v，  result 计算后的结果
void A_plus_vector(double v[SIZE],double result[SIZE])
{
    double a[SIZE];
    int i=0;
    for(i=0;i<SIZE;i++)
    {
        //a 的公式中，i从1开始，c语言是从0开始，故而公式上加了1. 另外sin(x),x是弧度
        a[i]=(1.64-0.024*(i+1))*sin(0.2*(i+1))-0.64*exp(0.1/(i+1));

    }

    double b=0.16;
    double c=-0.064;
    //A的第一行乘以v
    result[0]=a[0]*v[0]+b*v[1]+c*v[2];
    //A的第二行乘以v
    result[1]=b*v[0]+a[1]*v[1]+b*v[2]+c*v[3];
    //A的第三行到A的第499行乘v
    for(i=2;i<499;i++)
    {
        result[i]=v[i-2]*c+v[i-1]*b+v[i]*a[i]+v[i+1]*b+v[i+2]*c;
    }

    //A的第500行乘v
    result[499]=c*v[497]+b*v[498]+a[499]*v[499]+b*v[500];
    //A的第501行乘以v
    result[500]=c*v[498]+b*v[499]+a[500]*v[500];

}

//两个向量相乘，对应公式3
double dot_multiply(double v1[SIZE],double v2[SIZE])
{
    double result=0;
    int i=0;
    for(i=0;i<SIZE;i++)
        result+=v1[i]*v2[i];
    return result;
}

//幂法
double power_method()
{
    double u0[SIZE];
    double u1[SIZE];
    double y0[SIZE];
    double y1[SIZE];
    int i=0;
    //初始化u0
    for(i=0;i<SIZE;i++)
        u0[i]=0.1;

    //计算y0和u1

    vector_normalize(u0,y0);
    A_plus_vector(y0,u1);
    double beta1=dot_multiply(y0,u1);


    //计算y1和u2


    vector_normalize(u1,y1);

    A_plus_vector(y1,u0);
    double beta2=dot_multiply(y1,u0);



    //循环计算，直到达到误差精度。
    while((fabs(beta2-beta1)/fabs(beta2))>pow(10,-12))
    {

        vector_normalize(u0,y0);
        A_plus_vector(y0,u1);
        beta1=dot_multiply(y0,u1);

        vector_normalize(u1,y1);

        A_plus_vector(y1,u0);
        beta2=dot_multiply(y1,u0);

    };

    return beta2;


}




/*********************************************************
带原点位移反幂法部分

由于要求A*u=y,故使用Gauss消去法。同时可以求Det(A)
**********************************************************/

//对A进行p的位移，即对(A-pI)u=y*x进行高斯消去，返回u和Det(A)。
//p 位移量,y是公式中的y， result为计算后的结果，返回值为Det(A)
double Gauss_Elimination(double p,double y[SIZE],double result[SIZE])
{
    double a[501];
    int i=0;
    for(i=0;i<501;i++)
    {
        //a 的公式中，i从1开始，c语言是从0开始，故而公式上加了1. 另外sin(x),x是弧度
        a[i]=(1.64-0.024*(i+1))*sin(0.2*(i+1))-0.64*exp(0.1/(i+1))-p;//p是位移量
    }
    double b=0.16;
    double c=-0.064;


    /****************************************************
    由于矩阵的特殊性，每次只需要对下面三行做Gauss消去
    为了不存储0结果。定义一个501*5的矩阵。五列为
    c b ai  b c或者他们消去后的结果。
    ******************************************************/
    int j=0;
    double temp_matrix[SIZE][5];
    //对temp_matrix进行赋值

    for(i=0;i<SIZE;i++)
    {
        temp_matrix[i][0]=c;
        temp_matrix[i][1]=b;
        temp_matrix[i][2]=a[i];
        temp_matrix[i][3]=b;
        temp_matrix[i][4]=c;
        result[i]=y[i];
    }
    double m=0;

    //消去过程，只需要消去i下面两行
    for(i=0;i<499;i++)
    {
        //消去i下面第一行
        m=temp_matrix[i+1][1]/temp_matrix[i][2];
        for(j=1;j<4;j++)
        {
            temp_matrix[i+1][j]=temp_matrix[i+1][j]-m*temp_matrix[i][j+1];
        }

        result[i+1]=result[i+1]-m*result[i];

        //消去i下面第二行
        m=temp_matrix[i+2][0]/temp_matrix[i][2];
        for(j=0;j<3;j++)
        {
            temp_matrix[i+2][j]=temp_matrix[i+2][j]-m*temp_matrix[i][j+2];
        }

        result[i+2]=result[i+2]-m*result[i];
    }
    //消去最后一行
    m=temp_matrix[500][1]/temp_matrix[499][2];
    for(j=1;j<3;j++)
    {
        temp_matrix[500][j]=temp_matrix[500][j]-m*temp_matrix[499][j+1];
    }
    result[500]=result[500]-m*result[499];



    //回代过程，同样只需要回代前面两个
    for(i=500;i>1;i--)
    {
        m=temp_matrix[i-1][3]/temp_matrix[i][2];
        temp_matrix[i-1][3]=temp_matrix[i-1][3]-m*temp_matrix[i][2];
        result[i-1]=result[i-1]-m*result[i];

        m=temp_matrix[i-2][4]/temp_matrix[i][2];
        temp_matrix[i-2][4]=temp_matrix[i-2][4]-m*temp_matrix[i][2];
        result[i-2]=result[i-2]-m*result[i];
    }
    m=temp_matrix[0][3]/temp_matrix[1][2];
    temp_matrix[0][3]=temp_matrix[0][3]-m*temp_matrix[1][2];
    result[0]=result[0]-m*result[1];

    double Det=1;
    for(i=0;i<SIZE;i++)
    {
        Det=Det*temp_matrix[i][2];
        result[i]=result[i]/temp_matrix[i][2];
    }

    return Det;

}


//带位移的反幂法
double anti_power_method(double p)
{
    double u0[SIZE];
    double u1[SIZE];
    double y0[SIZE];
    double y1[SIZE];
    int i=0;
    //初始化u0
    for(i=0;i<SIZE;i++)
        u0[i]=1;

    //计算y0和u1

    vector_normalize(u0,y0);
    Gauss_Elimination(p,y0,u1);
    double beta1=dot_multiply(y0,u1);


    //计算y1和u2


    vector_normalize(u1,y1);
    Gauss_Elimination(p,y1,u0);
    double beta2=dot_multiply(y1,u0);



    //循环计算，直到达到误差精度。
    while((fabs(beta2-beta1)/fabs(beta2))>pow(10,-12))
    {

        vector_normalize(u0,y0);
        Gauss_Elimination(p,y0,u1);
        beta1=dot_multiply(y0,u1);

        vector_normalize(u1,y1);

        Gauss_Elimination(p,y1,u0);
        beta2=dot_multiply(y1,u0);

    };

    return 1/beta2+p;
}

int main()
{
    double lambda_max=power_method();
    double lambda1,lambda501;
    //如果lambda>0,则lambda501=lambda,且为按模最大。则lambda1为最接近-lambda501的
    printf("\n\n\n******************************\nanswer of question 1 is:\n");
    if(lambda_max>0)
    {
        lambda501=lambda_max;
        printf("lambda501=%.14e\n",lambda501);
        lambda1=anti_power_method(-lambda501);
        printf("lambda1=%.14e\n",lambda1);
    }


    //如果lambda<0,则lambda1=lambda，且为按摩最大。则lambda501为最接近-lambda1的
    if(lambda_max<0)
    {
        lambda1=lambda_max;
        printf("lambda1=%.14e\n",lambda1);
        lambda501=anti_power_method(-lambda1);
        printf("lambda501=%.14e\n",lambda501);
    }

    double lambdas=anti_power_method(0);
    printf("lambdas=%.14e\n",lambdas);

    //求于uk=lambda1+k(lambda501-lambda1)/40 最接近的lambda
    printf("\n\n\n******************************\nanswer of quetion 2 is:\n");
    int k=0;
    double uk=0;
    for(k=1;k<40;k++)
    {
        uk=lambda1+k*(lambda501-lambda1)/40;
        double lambda=anti_power_method(uk);
        printf("lambda_i%d=%.14e\n",k,lambda);
    }
    printf("\n\n\n******************************\nanswer of quetion 3 is:\n");
    //求条件数：对于非奇异实对称矩阵，cond为按模最大特征值/按模最小特征值的绝对值
    printf("cond(A)2=%.14e\n",fabs(lambda_max/lambdas));
    //求Det ，直接调用Gauss消去法
    double y[SIZE];
    double result[SIZE];
    double Det=Gauss_Elimination(0,y,result);
    printf("Det(A)=%.14e\n",Det);

    return 1;

}
