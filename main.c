#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void entering_func(int n, int mp[])
{
    int i;
    for (i = 0; i < n+1 ; i++)
    {
        printf("Please enter the %d. degree element\n", i);

        scanf("%d", &mp[i]);                                            //Belirtilen fonksiyonun katsayi degerlerini sakladigimiz dizi.

    }
}
void printing_func(int n, int mp[])
{
    int i;
    printf("Entered function is \n");

    for (i = n; i >= 0; i--) {
        if (i != 0) {
            printf("%d(x^%d) + ", mp[i], i);            //Katsayili bir sekilde fonksiyonun ciktisinin alinmasi.
        }
        else
            printf("%d", mp[i]);
    }
    printf("\n");
}
double func_resulter(double e,int mp[], int n)
{
    int i;
    double sum=0;

    for (i = 0; i <= n ; i++)
    {
        sum += pow(e,i) * mp[i];                    //Katsayilar dizisinden yararlanarak belirtilen deger icin girilen fonksiyonun sonucu hesaplayan fonksiyon.
    }

    return sum;
}
void entering_edges(double *f, double *l)
{
    double ff,ll;
    printf("Please enter the first and last edges of function.\n");
    printf("Please enter the first edge. \n");
    scanf("%lf", &ff);
    printf("Please enter the last edge. \n");
    scanf("%lf", &ll);                                                                  //Bazi islemler icin gerekli olan sinir degerlerinin alinamsi.

    *f = ff;
    *l = ll;
}


void bisection(double f, double l,int n, int mp[]){
    printf("\n\n");
    int iterationCount = 0;
    double sum1,sum2,error,epsilon = 99999,result;



    printf("What is the epslion ?\n");
    scanf("%lf", &epsilon);

    if(func_resulter(f, mp, n) == 0)
    {
        printf("Root is = %lf", f);
    }
    else if(func_resulter(l, mp, n) == 0)
    {
        printf("Root is = %lf", l);
    }
    else if (func_resulter((f + l) / 2,mp,n) == 0)
    {
        printf("\n\nIteration Has stopped, Root has been found. \nIteration Count = %d, Error = %lf, Root %lf\n\n\n", 1, 0.0, (f+l)/2);
    }

    else
    {
        do{
            iterationCount++;

            sum1 = func_resulter(f,mp,n);
            sum2 = func_resulter(l,mp,n);


            result = sum1 * sum2;
            error = ((l-f) / pow(2,iterationCount));

            if (result < 0 && error > epsilon){
                if (sum2 * func_resulter((f + l) / 2,mp,n) > 0) {
                    l = (f + l) / 2;
                }
                else {
                    f = (f + l) / 2;
                }
            }

        }while(error > epsilon);

        printf("\n\nIteration Has stopped, Root has been found. \nIteration Count = %d, Error = %lf, Root %lf\n\n\n", iterationCount, error, (f+l)/2);
        printf("\n\n");
    }







}
void regula_falsi(double f, double l,int n, int mp[]){
    int iterationCount = 1;
    double sum1,sum2,error,epsilon,edge;
    printf("What is the epslion ? \n");
    scanf("%lf", &epsilon);



    do{
        iterationCount++;


        sum1 = func_resulter(f,mp,n);
        sum2 = func_resulter(l,mp,n);

        edge = ((f * sum2) - (l * sum1)) / (sum2 - sum1);

        if (func_resulter(edge, mp, n) * func_resulter(f, mp, n)  < 0) {
            l = edge;
        }
        else {
            f = edge;
        }

        error = (l - f) / pow(2,iterationCount);
    }while(error > epsilon);

    printf("\n\nIteration Has stopped, Root has been found. \nIteration Count = %d, Error = %lf, Root %lf\n\n\n", iterationCount, error, edge);
    printf("\n\n");
}

void newtonRaphson(int n,int mp[]){
    int *d_mp;
    double val,epsilon,error,edge;

    printf("Please enter the begining value. \n");
    scanf("%lf", &val);

    printf("Please enter the epsilon\n");
    scanf("%lf", &epsilon);

    d_mp = (int*) malloc(n * sizeof(int));

    if (d_mp == NULL)
    {
        printf("Error on memorry allocation !\n");
        exit(0);
    }

    printf("Please enter the derivative of the function\n");
    entering_func(n-1,d_mp);
    printing_func(n-1,d_mp);


    do{
        edge = val - ((func_resulter(val,mp,n))/func_resulter(val,d_mp,n-1));
        error = fabs(edge - val);
        val = edge;

    }while(error > epsilon);
    printf("Root is approximately %lf", val);


}

void entering_matrix(double matrix[50][50], int row, int column){

    printf("Please enter the elements of matrix.\n");
    for (int i = 1; i <= row; i++)
    {
        for (int j = 1; j <= column; j++)
        {
            printf("Please enter the [%d][%d]\n",i,j);
            scanf("%lf", &matrix[i][j]);
        }
    }




}
void printing_matrix(double matrix[50][50], int row, int column){
    for (int i = 1; i<=row;i++){
        for(int j = 1; j <= column; j++)
        {
            printf("%.3lf \t", matrix[i][j]);
        }
        printf("\n");
    }
}

double determinant(double matrix[50][50], int dimension)
{
    /*
        -   Asagidaki kodun calisma mantigi elimizdeki matrisi recursive sekilde kucuk hallere bolerek kucuk parcalarin
        -   determinantini bulmaya yoneliktir.
        -   https://www.thecalculator.co/includes/forms/assets/img/Matrix%20determinant%204x4%20formula.jpg
        -   Islemin daha rahat yapilmasi adina matrisler 1. indisten baslamistir.

     */
    double s = 1, det = 0, temp[50][50];

    if(dimension == 1)
    {
        return matrix[1][1];
    }

    else
    {
        int k, l;
        for (int x = 1; x <= dimension; x++)
        {
            k = 1;
            l = 1;

            for (int i = 1; i <= dimension; i++)
            {
                for (int j = 1; j <= dimension; j++)
                {
                    temp[i][j] = 0;
                    if (i != 1 && j != x)
                    {
                        temp[k][l] = matrix[i][j];
                        if (l <= (dimension - 2))
                        {
                            l++;
                        }
                        else
                        {
                            l = 1;
                            k++;
                        }
                    }
                }
            }

            det = det + s * (matrix[1][x] * determinant(temp, k - 1));
            s *= -1;
        }
    }

    return det;

}
void transpose(double inverse[50][50], double matrix[50][50], int dimension)
{

    for (int i = 1; i <= dimension; i++)
    {
        for (int j = 1; j <= dimension; j++)
        {
            inverse[i][j] = matrix[j][i];
        }
    }

}
void cofactor(double inverse[50][50], double matrix[50][50], int dimension)
{
    /*
     * Kofaktor matristen yararlanarak bir islem yapiliyor.
     * https://wikimedia.org/api/rest_v1/media/math/render/svg/42fc09808599cf39611b7dce5fc27aa6e7bd424c
     * Buradaki verilere ulasmak icin de determinanttan yararlaniliyor.
    */

    double temp[50][50], cof[50][50];
    for (int k = 1; k <= dimension; k++)
    {
        for (int l = 1; l <= dimension; l++)
        {
            int m = 1;
            int n = 1;

            for (int i = 1; i <= dimension; i++)
            {
                for (int j = 1; j <= dimension; j++)
                {
                    if (i != k && j != l)
                    {
                        temp[m][n] = matrix[i][j];
                        if (n <= (dimension - 2))
                        {
                            n++;
                        }
                        else
                        {
                            n = 1;
                            m++;
                        }
                    }
                }
            }
            printf("\n");
            //printing_matrix(temp, dimension,dimension);printf("M = %d", m);
            //printf("\n");
            cof[k][l] = pow(-1, k + l) * determinant(temp, dimension - 1);
            //printf("DETT = %lf", determinant(temp, dimension - 1));
        }
    }


    transpose(inverse, cof, dimension);
}
void inverse_matrix(double matrix[50][50], int dimension)
{
    /*
     * Elimizdeki kofaktor matrisi matrisin determinantina boldugumuzde cikan sonuc bize matrisin
     * tersini verdiginden onceden hazirladigimiz cofactor dizisini determinanta boluyoruz.
     * https://www.onlinemath4all.com/images/inverseofamatrix.png
     * adj(A) = cof(A) ^ t
    */

    double inverse[50][50];

    cofactor(inverse, matrix, dimension);

    for (int i = 1; i <= dimension; i++)
    {
        for (int j = 1; j <= dimension; j++)
        {
            inverse[i][j] /= determinant(matrix,dimension);
        }
    }

    printf("Inverted matrix is down below!\n\n");
    printing_matrix(inverse, dimension, dimension);
}

void gauss_elemination(double matrix[50][50], int n)
{
    /*
     * Alt ucgen matris elde ediyoruz.
     * https://www.youtube.com/watch?v=jOC4fMgl7TI
     * Ettigimiz alt ucgeni sonra birim ucgene cevirip en son degerleri tek boyutlu dizinin icine atip
     * ciktisini aliyoruz.
    */

    double ratio , temp[50];
    for (int i = 1; i <= n ; i++)
    {
        for(int j = i+1; j <= n; j++ )
        {
            ratio = matrix[j][i] / matrix[i][i];
            for(int k = 1; k <= n+1 ; k++)
            {
                matrix[j][k] = matrix[j][k] - ratio * matrix[i][k];
            }
        }
    }

    temp[n] = matrix[n][n+1]/matrix[n][n];

    for(int i = n-1; i >= 1; i--)
    {
        temp[i] = matrix[i][n+1];
        for(int j = i+1;j <= n; j++)
        {
            temp[i] = temp[i] - matrix[i][j] * temp[j];
        }
    temp[i] = temp[i]/matrix[i][i];
    }


    for(int i =1; i <= n; i++)
    {
        printf("x[%d] = %0.3f\n",i, temp[i]);
    }


}


void gs_fx(double matrix[50][50],int n,double mp[],int x,double deltas[])
{
    double sum;
    sum = matrix[x][n+1];
    for (int i = 1; i <= n ; i++)
    {
        if (n+1-i != x)
        {
            sum -= matrix[x][n+1 - i] * mp[n+1-i];
        }
    }
    sum /= matrix[x][x];
    deltas[x] = fabs(sum - mp[x]);
    mp[x] = sum;

}
void gauss_seidel(double matrix[50][50],int n)
{

    int x = 1, flag = 0;
    double mp[50], deltas[50],epsilon;

    printf("What is the value of epslion ? \n");
    scanf("%lf", &epsilon);

    for (int i = 1; i <= n; ++i)
    {
        deltas[i]= 9999;
        printf("Please enter the X%d. variable's starting value. \n",i);
        scanf("%lf", &mp[i]);
    }

    while(flag == 0)
    {

        gs_fx(matrix, n, mp,x,deltas);
        x++;
        if(x > n){
            x = 1;
        }

        flag = 1;


        int y = 1;
        while(y <= n && flag == 1)
        {
            if(deltas[y] > epsilon)
            {
                flag = 0;
            }
            y++;
        }
    }

    for (int i = 1; i <= n; i++)
    {
        printf("X%d = %lf\n", i, round(mp[i]));
    }
}

double forward_det(int n, int mp[], double h, double x)
{
    double sum;
    sum = func_resulter(x + h, mp, n) - func_resulter(x, mp, n);
    sum /= h;
    return sum;
}
double backward_det(int n, int mp[], double h, double x)
{
    double sum;
    sum = func_resulter(x, mp, n) - func_resulter(x - h, mp, n);
    sum /= h;
    return sum;
}
double centeral_det(int n, int mp[], double h, double x)
{
    double sum;
    sum = func_resulter(x + h, mp, n) - func_resulter(x - h, mp, n);
    sum /= 2*h;
    return sum;
}
double detirative(int n, int mp[])
{
    int theV;
    double h,x;

    printf("Which detirative operation do you want to do ?\n");
    printf("1 - Forward detirative \n2 - Backward detirative\n3 - Centeral detirative\n");
    scanf("%d", &theV);
    printf("What is the value of h ?\n");
    scanf("%lf", &h);
    printf("What is the value of x ?\n");
    scanf("%lf", &x);

    if (theV == 1)
    {
        return forward_det(n, mp, h, x);
    }
    else if (theV == 2)
    {
        return backward_det(n, mp, h, x);
    }
    else if(theV == 3)
    {
        return centeral_det(n, mp, h, x);
    }
    else
        exit(2);
}

double sum_trapez(int n, int mp[], double x0,double dx, int x)
{
    double sum = 0;

    for (int i = 1; i < x- 1 ; i++ )
    {
        sum += func_resulter(x0 + i*dx, mp, n);
    }

    return sum;
}
double trapez(int n, int mp[])
{
    int x;
    double sum, x0,xn, dx;

    printf("\nPlease enter the value 'n'.\n");
    scanf("%d",&x);

    printf("Please enter the value x0. \n");
    scanf("%lf", &x0);

    printf("Please enter the value xn \n");
    scanf("%lf", &xn);

    dx = ( xn - x0 ) / x;
    sum = dx * ((( func_resulter(x0, mp, n) + func_resulter(xn, mp, n) ) / 2 ) + sum_trapez(n, mp, x0, dx, x));

    return sum;



}

double sum_Simpson13(int n, int mp[], double x0,double dx, int x, int decider)
{
    double sum = 0;

    if (decider == 0)
    {
        for (int i = 1; i < x-1; i = i + 2)
        {
            sum += func_resulter(x0 + i*dx, mp, n);
        }
    }
    else
    {
        for (int i = 2; i < x-2; i = i + 2)
        {
            sum += func_resulter(x0 + i*dx, mp, n);
        }
    }

    return sum;
}
double simpson13(int n, int mp[])
{
    int x;
    double sum, x0,xn, dx;

    printf("\nPlease enter the value 'n'.(Remember! You have to enter an even number.)\n");
    scanf("%d",&x);

    printf("Please enter the value x0. \n");
    scanf("%lf", &x0);

    printf("Please enter the value xN. \n");
    scanf("%lf", &xn);

    dx = ( xn - x0 ) / x;

    sum = dx/3 * (func_resulter(x0, mp, n) + func_resulter(xn, mp, n) + 4 * sum_Simpson13(n, mp, x0, dx, x, 0) + 2 * sum_Simpson13(n, mp, x0, dx, x, 1));

    return sum;


}


void entering_gregory(int n, double x[], double y[50][50])
{

    for (int i = 0; i < n; i++)
    {
        printf("Please enter the value of x[%d] \n", i);
        scanf("%lf", &x[i]);
        printf("Please enter the value of y[%d]", i);
        scanf("%lf", &y[i][0]);
    }

}
void determining_dx(int n, double y[50][50], int *degree)
{
    int diff, i = 1;

    while (i < n && diff != 0)
    {
        for (int j = 0; j < n - i; j++)
        {
            y[j][i] = y [j + 1][i - 1] - y [j][i - 1];
        }
        for (int j = 0; j < n - i; j++)
        {
            if ( y[j][i - 1] == y[j + 1][i - 1] )
            {
                diff = 0;
            }
            else
            {
                diff = 1;
            }
        }
        i++;

    }

    *degree = i - 2;


    printf(" Matrix of dFx \n");

    for(i = 0; i < n; i++)
    {
        for(int j = 0; j < n-i ; j++)
        {
            printf("\t%f", y[i][j]);
        }
        printf("\n");
    }


}

int factorial ( int n )
{
    if ( n == 0 )
    {
        return 1;
    }
    else
    {
        return n * factorial(n-1);
    }
}
double main_gregory(int degree, double x[], double y[50][50], double a)
{

    double sum, result = 0;
    int j = 0;
    double h = x[1] - x [0];



    result = y[0][0];

    while (j < degree)
    {
        sum = 1;
        for(int i = 0; i <= j; i++)
        {
            sum *= (a - x[i]);
        }

        sum *= y[0][j + 1] / pow(h, j+1) / factorial(j + 1);
        j++;
        result += sum;

    }
    return result;

}
void gregory_newton()
{
    int n, degree;
    double result, a, x[50], y[50][50];


    printf("Please enter the number of tests. \n");
    scanf("%d", &n);

    printf("Please enter the value x \n");
    scanf("%lf", &a);


    entering_gregory(n,x,y);
    determining_dx(n,y,&degree);
    result = main_gregory(degree, x, y, a);

    printf("The result is = %lf", result);

}

int main() {

    int theV1 = 0,theV2 = 0;
    printf("Welcome to the numeric analysis calculator! \n");

    while(theV1 != -1){
        int *mp, n;

        printf("\n\n\n");
        printf("What do you want to do?\n");
        printf("Dou you want to do a functional problem(1) or matrixal operation(2)?\n");
        printf("For exiting the program, please enter value -1\n\n\n");
        scanf("%d",&theV1);

        if(theV1 != -1){
            if(theV1 == 1)
            {
                printf("Which opreation do you want to do?\n");
                printf("1 - Bisection \n2 - Regula-Falsi \n3 - Newton-Raphson \n4 - Detirative \n5 - Trapez \n6 - Simpson(1/3) \n7 - Gregory - Newton \n");
                scanf("%d", &theV2);

                if(theV2 != 7 )
                {
                    printf("What is the degree of function? \n");
                    scanf("%d", &n);
                    mp = (int*) malloc(n+1 * sizeof(int));

                    if (mp == NULL)
                    {
                        printf("Error on memorry allocation !\n");
                        exit(0);
                    }
                    entering_func(n,mp);
                    printing_func(n,mp);
                }

                if (theV2 == 1)
                {
                    double f,l;
                    entering_edges(&f,&l);
                    bisection(f,l,n,mp);
                }
                else if (theV2 == 2)
                {
                    double f,l;
                    entering_edges(&f,&l);
                    regula_falsi(f,l,n,mp);

                }
                else if (theV2 == 3)
                {
                    newtonRaphson(n,mp);
                }
                else if(theV2 == 4)
                {
                    printf("The detirative is = %lf\n",detirative(n,mp));
                }
                else if(theV2 == 5)
                {
                    printf("Integral of your func in trapez is = %lf\n", trapez(n,mp));
                }
                else if(theV2 == 6)
                {
                    printf("Integral of your func in simpson(1/3) is %lf\n", simpson13(n,mp));
                }
                else if (theV2 == 7)
                {
                    gregory_newton();
                }
            }
            else if (theV1 == 2){
                double matrix[50][50];
                int row,column;


                printf("Which opreation do you want to do?\n");
                printf("1- Inversing a matrix\n2- Gauss Elemination\n3- Gauss Seidal");
                scanf("%d", &theV2);

                printf("Please enter your matrix's dimensions. Firstly column then row.\n");
                scanf("%d",&column);
                scanf("%d",&row);

                if(theV2 == 1)
                {
                    entering_matrix(matrix, row, column);
                    printing_matrix(matrix,row, column);
                    //printf("Det = %lf",determinant(matrix,row));
                    inverse_matrix(matrix,row);
                }
                else if (theV2 == 2)
                {
                    entering_matrix(matrix,column,row+1);
                    printing_matrix(matrix,column,row+1);
                    gauss_elemination(matrix,row);
                }
                else if(theV2 == 3)
                {
                    printf("Please enter the matrix in diagonally dominant form \n");
                    entering_matrix(matrix, column, row+1);
                    printing_matrix(matrix,column,row+1);
                    gauss_seidel(matrix, row);
                }
            }
        }
    }

    return 0;
}

