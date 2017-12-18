#include <stdio.h>
#include <stdlib.h>
#include<math.h>

#define SIZE 501

//ע����������Ŀ�и�����sin(0.2i)�У�û��ע���ǻ��Ȼ��ǽǶȣ����㷨ȫ�����ջ��ȼ��㡣


/***************************************************************
�ݷ�����
*************************************************************/

// �����������ڵ�ģ������Ӧ���еĵ�һ���͵ڶ�����ʽ,��������һ��
//v ����һ��������   nomalized_v  ��һ��֮��Ľ��
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

//A�����������v,����Ҫ��A�ķ�0���ݶ������棬�ʵ����������ֳ˷������з�0Ԫ�ؾ�������
//v  A*v�е�v��  result �����Ľ��
void A_plus_vector(double v[SIZE],double result[SIZE])
{
    double a[SIZE];
    int i=0;
    for(i=0;i<SIZE;i++)
    {
        //a �Ĺ�ʽ�У�i��1��ʼ��c�����Ǵ�0��ʼ���ʶ���ʽ�ϼ���1. ����sin(x),x�ǻ���
        a[i]=(1.64-0.024*(i+1))*sin(0.2*(i+1))-0.64*exp(0.1/(i+1));

    }

    double b=0.16;
    double c=-0.064;
    //A�ĵ�һ�г���v
    result[0]=a[0]*v[0]+b*v[1]+c*v[2];
    //A�ĵڶ��г���v
    result[1]=b*v[0]+a[1]*v[1]+b*v[2]+c*v[3];
    //A�ĵ����е�A�ĵ�499�г�v
    for(i=2;i<499;i++)
    {
        result[i]=v[i-2]*c+v[i-1]*b+v[i]*a[i]+v[i+1]*b+v[i+2]*c;
    }

    //A�ĵ�500�г�v
    result[499]=c*v[497]+b*v[498]+a[499]*v[499]+b*v[500];
    //A�ĵ�501�г���v
    result[500]=c*v[498]+b*v[499]+a[500]*v[500];

}

//����������ˣ���Ӧ��ʽ3
double dot_multiply(double v1[SIZE],double v2[SIZE])
{
    double result=0;
    int i=0;
    for(i=0;i<SIZE;i++)
        result+=v1[i]*v2[i];
    return result;
}

//�ݷ�
double power_method()
{
    double u0[SIZE];
    double u1[SIZE];
    double y0[SIZE];
    double y1[SIZE];
    int i=0;
    //��ʼ��u0
    for(i=0;i<SIZE;i++)
        u0[i]=0.1;

    //����y0��u1

    vector_normalize(u0,y0);
    A_plus_vector(y0,u1);
    double beta1=dot_multiply(y0,u1);


    //����y1��u2


    vector_normalize(u1,y1);

    A_plus_vector(y1,u0);
    double beta2=dot_multiply(y1,u0);



    //ѭ�����㣬ֱ���ﵽ���ȡ�
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
��ԭ��λ�Ʒ��ݷ�����

����Ҫ��A*u=y,��ʹ��Gauss��ȥ����ͬʱ������Det(A)
**********************************************************/

//��A����p��λ�ƣ�����(A-pI)u=y*x���и�˹��ȥ������u��Det(A)��
//p λ����,y�ǹ�ʽ�е�y�� resultΪ�����Ľ��������ֵΪDet(A)
double Gauss_Elimination(double p,double y[SIZE],double result[SIZE])
{
    double a[501];
    int i=0;
    for(i=0;i<501;i++)
    {
        //a �Ĺ�ʽ�У�i��1��ʼ��c�����Ǵ�0��ʼ���ʶ���ʽ�ϼ���1. ����sin(x),x�ǻ���
        a[i]=(1.64-0.024*(i+1))*sin(0.2*(i+1))-0.64*exp(0.1/(i+1))-p;//p��λ����
    }
    double b=0.16;
    double c=-0.064;


    /****************************************************
    ���ھ���������ԣ�ÿ��ֻ��Ҫ������������Gauss��ȥ
    Ϊ�˲��洢0���������һ��501*5�ľ�������Ϊ
    c b ai  b c����������ȥ��Ľ����
    ******************************************************/
    int j=0;
    double temp_matrix[SIZE][5];
    //��temp_matrix���и�ֵ

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

    //��ȥ���̣�ֻ��Ҫ��ȥi��������
    for(i=0;i<499;i++)
    {
        //��ȥi�����һ��
        m=temp_matrix[i+1][1]/temp_matrix[i][2];
        for(j=1;j<4;j++)
        {
            temp_matrix[i+1][j]=temp_matrix[i+1][j]-m*temp_matrix[i][j+1];
        }

        result[i+1]=result[i+1]-m*result[i];

        //��ȥi����ڶ���
        m=temp_matrix[i+2][0]/temp_matrix[i][2];
        for(j=0;j<3;j++)
        {
            temp_matrix[i+2][j]=temp_matrix[i+2][j]-m*temp_matrix[i][j+2];
        }

        result[i+2]=result[i+2]-m*result[i];
    }
    //��ȥ���һ��
    m=temp_matrix[500][1]/temp_matrix[499][2];
    for(j=1;j<3;j++)
    {
        temp_matrix[500][j]=temp_matrix[500][j]-m*temp_matrix[499][j+1];
    }
    result[500]=result[500]-m*result[499];



    //�ش����̣�ͬ��ֻ��Ҫ�ش�ǰ������
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


//��λ�Ƶķ��ݷ�
double anti_power_method(double p)
{
    double u0[SIZE];
    double u1[SIZE];
    double y0[SIZE];
    double y1[SIZE];
    int i=0;
    //��ʼ��u0
    for(i=0;i<SIZE;i++)
        u0[i]=1;

    //����y0��u1

    vector_normalize(u0,y0);
    Gauss_Elimination(p,y0,u1);
    double beta1=dot_multiply(y0,u1);


    //����y1��u2


    vector_normalize(u1,y1);
    Gauss_Elimination(p,y1,u0);
    double beta2=dot_multiply(y1,u0);



    //ѭ�����㣬ֱ���ﵽ���ȡ�
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
    //���lambda>0,��lambda501=lambda,��Ϊ��ģ�����lambda1Ϊ��ӽ�-lambda501��
    printf("\n\n\n******************************\nanswer of question 1 is:\n");
    if(lambda_max>0)
    {
        lambda501=lambda_max;
        printf("lambda501=%.14e\n",lambda501);
        lambda1=anti_power_method(-lambda501);
        printf("lambda1=%.14e\n",lambda1);
    }


    //���lambda<0,��lambda1=lambda����Ϊ��Ħ�����lambda501Ϊ��ӽ�-lambda1��
    if(lambda_max<0)
    {
        lambda1=lambda_max;
        printf("lambda1=%.14e\n",lambda1);
        lambda501=anti_power_method(-lambda1);
        printf("lambda501=%.14e\n",lambda501);
    }

    double lambdas=anti_power_method(0);
    printf("lambdas=%.14e\n",lambdas);

    //����uk=lambda1+k(lambda501-lambda1)/40 ��ӽ���lambda
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
    //�������������ڷ�����ʵ�Գƾ���condΪ��ģ�������ֵ/��ģ��С����ֵ�ľ���ֵ
    printf("cond(A)2=%.14e\n",fabs(lambda_max/lambdas));
    //��Det ��ֱ�ӵ���Gauss��ȥ��
    double y[SIZE];
    double result[SIZE];
    double Det=Gauss_Elimination(0,y,result);
    printf("Det(A)=%.14e\n",Det);

    return 1;

}
