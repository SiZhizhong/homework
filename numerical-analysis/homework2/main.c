#include <stdio.h>
#include <stdlib.h>
#include<math.h>

#define SIZE 10

//符号函数
double sgn(double s)
{
    if(s<0)
        return -1;
    else
        return 1;
}
//矩阵乘以向量，输出矩阵
void matrix_multipy_vector(double input_matrix[SIZE][SIZE],double input_vector[SIZE],double output_vector[SIZE])
{
    int i,j;
    for(i=0;i<SIZE;i++)
    {
        output_vector[i]=0;
        for(j=0;j<SIZE;j++)
            output_vector[i]=output_vector[i]+input_matrix[i][j]*input_vector[j];
    }
}

//向量点积，输出标量
double vector_dot_multiply_vector(double input_vector1[SIZE],double input_vector2[SIZE])
{
    double result=0;
    int i=0;
    for(i=0;i<SIZE;i++)
        result+=input_vector1[i]*input_vector2[i];
    return result;
}

//纵向量乘横向量，返回矩阵
void vector_multiply(double input_vector1[SIZE],double input_vector2[SIZE],double output_matrix[SIZE][SIZE])
{
    int i=0,j=0;
    for(i=0;i<SIZE;i++)
    {
        for(j=0;j<SIZE;j++)
            output_matrix[i][j]=input_vector1[i]*input_vector2[j];
    }
}

//向量减法，A-t*B返回向量
void vector_sub(double input_vector1[SIZE],double input_vector2[SIZE],double t,double output_vector[SIZE])
{
    int i;
    for(i=0;i<SIZE;i++)
        output_vector[i]=input_vector1[i]-t*input_vector2[i];
}

//矩阵减法,返回矩阵
void matrix_sub(double input_matrix1[SIZE][SIZE],double input_matrix2[SIZE][SIZE],double output_matrix[SIZE][SIZE])
{
    int i=0,j=0;
    for(i=0;i<SIZE;i++)
    {
        for(j=0;j<SIZE;j++)
            output_matrix[i][j]=input_matrix1[i][j]-input_matrix2[i][j];
    }
}

//矩阵转置
void matrix_transpose(double input_matrix[SIZE][SIZE],double output_matrix[SIZE][SIZE])
{
    int i,j;
    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
        output_matrix[i][j]=input_matrix[j][i];
}

//向量除标量
void vector_div_scalar(double input_vector[SIZE],double s,double output_vector[SIZE])
{
    int i;
    for(i=0;i<SIZE;i++)
        output_vector[i]=input_vector[i]/s;
}


//矩阵乘法
void matrix_multiply(double input_matrix1[SIZE][SIZE],double input_matrix2[SIZE][SIZE],double output_matrix[SIZE][SIZE])
{
    int i,j,k;
    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
        {
            output_matrix[i][j]=0;
            for(k=0;k<SIZE;k++)
                output_matrix[i][j]+=input_matrix1[i][k]*input_matrix2[k][j];
        }
}


//对输入矩阵进行拟上三角化
void Hessenberg( double input_matrix[SIZE][SIZE], double output_matrix[SIZE][SIZE])
{
    int i,j,k,all_zero_flag;
    double dr,cr,hr,tr;
    double ur[SIZE], pr[SIZE],qr[SIZE],wr[SIZE];
    double temp_matrix[SIZE][SIZE];
    double temp_matrix_trans[SIZE][SIZE];
    double v_multi_v_matirx[SIZE][SIZE];
    //利用中间矩阵进行计算，防止破坏原来输入矩阵；
    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
        temp_matrix[i][j]=input_matrix[i][j];


    //对矩阵的每一列做运算
    for (j=0;j<SIZE-2;j++)
    {
        //判断是不是所有元素都是0
        for (i=j+2;i<SIZE;i++)
           all_zero_flag=all_zero_flag*(temp_matrix[i][j]==0);
        if(all_zero_flag==1) continue;//如果全部为0，跳到书中P62页步骤5
        else
        {
            dr=0;
            for (i=j+1;i<SIZE;i++)
                dr+=temp_matrix[i][j]*temp_matrix[i][j];
            dr=sqrt(dr);
            if(temp_matrix[j+1][j]==0)
                cr=dr;
            else
                cr=-sgn(temp_matrix[j+1][j])*dr;
            hr=cr*cr-cr*temp_matrix[j+1][j];
            for(k=0;k<j+1;k++)
                ur[k]=0;
            ur[j+1]=temp_matrix[j+1][j]-cr;
            for(k=j+2;k<SIZE;k++)
                ur[k]=temp_matrix[k][j];
            matrix_transpose(temp_matrix,temp_matrix_trans);
            matrix_multipy_vector(temp_matrix_trans,ur,pr);
            vector_div_scalar(pr,hr,pr);
            matrix_multipy_vector(temp_matrix,ur,qr);
            vector_div_scalar(qr,hr,qr);
            tr=vector_dot_multiply_vector(pr,ur)/hr;
            vector_sub(qr,ur,tr,wr);
            vector_multiply(wr,ur,v_multi_v_matirx);
            matrix_sub(temp_matrix,v_multi_v_matirx,temp_matrix);
            vector_multiply(ur,pr,v_multi_v_matirx);
            matrix_sub(temp_matrix,v_multi_v_matirx,temp_matrix);

        }

    }

    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
        {
            output_matrix[i][j]=temp_matrix[i][j];
            if(fabs(output_matrix[i][j])<pow(10,-15))
                output_matrix[i][j]=0;
        }


}


void QR_decomposition(double input_matrix[SIZE][SIZE],double Q[SIZE][SIZE],double R[SIZE][SIZE])
{
    int i,j,k,all_zero_flag;
    double dr,cr,hr;
    double ur[SIZE], pr[SIZE],wr[SIZE];
    double temp_matrix[SIZE][SIZE];
    double temp_matrix_trans[SIZE][SIZE];
    double v_multi_v_matirx[SIZE][SIZE];
    //Q=I
    for (i=0;i<SIZE;i++)
    {
        for(j=0;j<SIZE;j++)
        {
            if(i==j) Q[i][j]=1;
            else Q[i][j]=0;
        }
    }

    //利用中间矩阵进行计算，防止破坏原来输入矩阵；
    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
        temp_matrix[i][j]=input_matrix[i][j];


    //对矩阵的每一列做运算
    for (j=0;j<SIZE-1;j++)
    {
        //判断是不是所有元素都是0
        for (i=j+1;i<SIZE;i++)
           all_zero_flag=all_zero_flag*(temp_matrix[i][j]==0);
        if(all_zero_flag==1) continue;//如果全部为0，跳到书中P62页步骤5
        else
        {
            dr=0;
            for (i=j;i<SIZE;i++)
                dr+=temp_matrix[i][j]*temp_matrix[i][j];
            dr=sqrt(dr);
            if(temp_matrix[j][j]==0)
                cr=dr;
            else
                cr=-sgn(temp_matrix[j][j])*dr;
            hr=cr*cr-cr*temp_matrix[j][j];
            for(k=0;k<j;k++)
                ur[k]=0;
            ur[j]=temp_matrix[j][j]-cr;
            for(k=j+1;k<SIZE;k++)
                ur[k]=temp_matrix[k][j];
            matrix_multipy_vector(Q,ur,wr);
            vector_div_scalar(ur,hr,ur);
            vector_multiply(wr,ur,v_multi_v_matirx);
            matrix_sub(Q,v_multi_v_matirx,Q);
            matrix_transpose(temp_matrix,temp_matrix_trans);
            matrix_multipy_vector(temp_matrix_trans,ur,pr);
            vector_div_scalar(ur,1/hr,ur);
            vector_multiply(ur,pr,v_multi_v_matirx);
            matrix_sub(temp_matrix,v_multi_v_matirx,temp_matrix);



        }

    }
    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
        {
            R[i][j]=temp_matrix[i][j];
            if(fabs(R[i][j])<pow(10,-15))
                R[i][j]=0;
        }


}

//用QR分解法，循环迭代
void QR_method(double input_matrix[SIZE][SIZE],double output_matrix[SIZE][SIZE])
{
    double R[SIZE][SIZE],Q[SIZE][SIZE];
    double Ak[SIZE][SIZE];
    int i,j;
    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
        Ak[i][j]=input_matrix[i][j];
    while(fabs(Ak[7][6])>pow(10,-20))
    {
        QR_decomposition(Ak,Q,R);
        matrix_multiply(R,Q,Ak);
    }
    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
        output_matrix[i][j]=Ak[i][j];

}

//二阶子阵，求特征值
void two_order_matrix(double a,double b,double c,double d,double *real,double *image)
{
    double p1,p2,p3,delta;
    p1=1;
    p2=-(a+d);
    p3=a*d-b*c;
    delta=p2*p2-4*p1*p3;
    if(delta>0) printf("error,2-order matrix give a real number eign value\n");
    else
    {
        *real=-p2/(2*p1);
        *image=sqrt(-delta)/(2*p1);
    }
}
//eign value 用二维数组保存，分别保存实部和虚部，如果是实数，这虚部为0
void eign_value(double input_vector[SIZE][SIZE],double eignvalue[SIZE][2])
{
    int i=0;
    double real,image;
    //不能算到0，因为不存在input_vector[0][0-1]这个单元
    for(i=SIZE-1;i>0;i--)
    {
        //得到一个实数特征值
        if(fabs(input_vector[i][i-1])<pow(10,-12))
        {
            eignvalue[i][0]=input_vector[i][i];
            eignvalue[i][1]=0;
        }
        else
        {
            //得到一对复数特征值
            two_order_matrix(input_vector[i-1][i-1],input_vector[i-1][i],input_vector[i][i-1],input_vector[i][i],&real,&image);
            eignvalue[i][0]=real;
            eignvalue[i][1]=image;
            eignvalue[i-1][0]=real;
            eignvalue[i-1][1]=-image;
            i--;
        }
    }
    if(i==0)
    {
        eignvalue[0][0]=input_vector[0][0];
        eignvalue[0][1]=0;
    }
}

//用列主要元素高斯消去法，在已知特征值的情况下，求特征向量
void eign_vector(double input_matrix[SIZE][SIZE],double lambda,double eignvector[SIZE])
{
    double A[SIZE][SIZE];
    int i,j,k,column_pivot_row_number;
    double temp,mik;
    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
        {
            if(i==j)
                 A[i][j]=input_matrix[i][j]-lambda;
            else
                 A[i][j]=input_matrix[i][j];
        }
    //消去过程
    for(i=0;i<SIZE-1;i++)
    {
        column_pivot_row_number=i;
        //选出列主元
        for(k=i+1;k<SIZE;k++)
            if(fabs(A[k][i])>fabs(A[column_pivot_row_number][i]))
            column_pivot_row_number=k;
        //交换
        for(j=0;j<SIZE;j++)
        {
            temp=A[i][j];
            A[i][j]=A[column_pivot_row_number][j];
            A[column_pivot_row_number][j]=temp;
        }
        for(k=i+1;k<SIZE;k++)
        {
            mik=A[k][i]/A[i][i];
            for(j=i;j<SIZE;j++)
                A[k][j]=A[k][j]-mik*A[i][j];
        }



    }
    //回代过程
        eignvector[SIZE-1]=1;
        for(i=SIZE-2;i>=0;i--)
        {
            temp=0;
            for(j=i+1;j<SIZE;j++)
                temp+=eignvector[j]*A[i][j];
            eignvector[i]=-temp/A[i][i];
        }

        //单位化特征向量
        temp=0;
        for(i=0;i<SIZE;i++)
            temp+=eignvector[i]*eignvector[i];
        temp=sqrt(temp);
        for(i=0;i<SIZE;i++)
            eignvector[i]=eignvector[i]/temp;




}

//数据太长，将数据保存到一个txt文档中
void matrix_print(double input_matrix[SIZE][SIZE],FILE *fp)
{

    int i=0,j=0;
    for (i=0;i<SIZE;i++)
    {
        fprintf(fp,"\nthe %dth row is:\n",i+1);
        for(j=0;j<SIZE;j++)
            fprintf(fp,"%.14e\t",input_matrix[i][j]);
        fprintf(fp,"\n");
    }
}

int main()
{
    double A[SIZE][SIZE],Hessenberg_A[SIZE][SIZE],Q[SIZE][SIZE],R[SIZE][SIZE];
    double QR_method_A[SIZE][SIZE];
    double eignvalue[SIZE][2];
    double eignvector[SIZE];

    int i,j;
    FILE *fp;
    fp=fopen("matrix.txt","a+");
    for(i=0;i<SIZE;i++)
        for(j=0;j<SIZE;j++)
        {
            if(i==j)
                A[i][j]=1.52*cos(i+1+1.2*j+1.2);//i+1,j+1
            else
                A[i][j]=sin(0.5*i+0.2*j+0.7);
        }
    fprintf(fp,"\n\n原矩阵A是:\n");
    matrix_print(A,fp);
    Hessenberg(A,Hessenberg_A);
    fprintf(fp,"\n\n拟上三角化矩阵A(n-1)是:\n");
    matrix_print(Hessenberg_A,fp);
    QR_decomposition(Hessenberg_A,Q,R);
    fprintf(fp,"\n\nA(n-1)QR分解Q是:\n");
    matrix_print(Q,fp);
    fprintf(fp,"\n\nA(n-1)QR分解R是:\n");
    matrix_print(R,fp);

    QR_method(Hessenberg_A,QR_method_A);
    fprintf(fp,"\n\nA(n-1)的QR迭代结束后矩阵是is:\n");
    matrix_print(QR_method_A,fp);

    eign_value(QR_method_A,eignvalue);
    fprintf(fp,"\n\n特征值是:\n");
    for(i=0;i<SIZE;i++)
    {
        if(eignvalue[i][1]==0)
            fprintf(fp,"特征值%.2d是\t%.14e\n",i+1,eignvalue[i][0]);
        else
        {
            fprintf(fp,"特征值%.2d是\t(%.14e) + (%.14e)i\n",i,eignvalue[i][0],eignvalue[i][1]);
        }
    }
    fprintf(fp,"\n\n实数特征值的特征向量是:\n\n");
    for(i=0;i<SIZE;i++)
    {
        //对实数根求特征向量
        if(eignvalue[i][1]==0)
        {

            eign_vector(A,eignvalue[i][0],eignvector);
            fprintf(fp,"特征值%.14e的特征向量是:\n",eignvalue[i][0]);

            for(j=0;j<SIZE;j++)
                fprintf(fp,"%.14e\t",eignvector[j]);
            fprintf(fp,"\n\n");

        }
    }
    return 0;
}
