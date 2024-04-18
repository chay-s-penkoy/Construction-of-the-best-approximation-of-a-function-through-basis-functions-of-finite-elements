#define _CRT_SECURE_NO_WARNINGS // для fopen 
#include <iostream>

#include <random>


#include <Dense>

#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;


using namespace std;

double my_func(double x) { // функция для приближения
    // sin(pow(x, 2)) - x * cos(x / 20);
    // 1 / (1 + pow(x, 2));
    return  sin(pow(x, 2)) - x * cos(x / 20);
}


void pogr(const int M, double* vector_c, double* x_metod, const char* name, FILE* errors) {

    double err_l_1 = 0.0; // асбсоллютные ошибки
    double err_l_2 = 0.0;
    double err_l_inf = 0.0;

    double val_l_1 = 0.0; // оценка значений f(x)
    double val_l_2 = 0.0;
    double val_l_inf = 0.0;



    for (int w = 0; w < M; w++) { // проходимся по всем интервалам

        double summ = vector_c[w];

        double chisl = x_metod[w];

        // абсоллютная погрешность
        err_l_1 += abs(summ - chisl);  // l1 
        err_l_2 += pow(summ - chisl, 2);  // l2
        if (abs(summ - chisl) > err_l_inf) { // l_inf
            err_l_inf = abs(summ - chisl);
        }

        // считаем знаменатели
        val_l_1 += abs(chisl);
        val_l_2 += pow(chisl, 2);
        if (abs(chisl) > val_l_inf) { // l_inf
            val_l_inf = abs(chisl);
        }

    }


    err_l_2 = sqrt(err_l_2);
    val_l_2 = sqrt(val_l_2);




    fprintf(errors, "\n");

    fprintf(errors, "|                         ВЫЧИСЛЕНИЕ ПОГРЕШНОСТЕЙ ДЛЯ %s                          | \n", name);
    fprintf(errors, "|               |       ||*||_1          |          ||*||_2         |        ||*||_inf         | \n");
    fprintf(errors, "|  абсоллютная  |    %14.8e      |     %14.8e       |     %14.8e       | \n", err_l_1, err_l_2, err_l_inf);
    fprintf(errors, "| относительная |    %14.8e      |     %14.8e       |     %14.8e       | \n", err_l_1 / val_l_1, err_l_2 / val_l_2, err_l_inf / val_l_inf);

 
}

void nevyazka(const int M, const int N, double* vector_b, double* matrix_A , double* x_metod, const char* name, FILE* errors) {

    double* nev; // 
    nev = new double[M];

    //double *A = (double*)matrix_A;

    //nev = new double [M];
    // матрично-векторное умножение 
    for (int i = 0; i < M; i++) {
        nev[i] = vector_b[i];
        for (int j = 0; (j < N) && (i + j < M); j++) {

            nev[i] -= matrix_A[ i * N + j] * x_metod[i + j];
            //Ax[i] += matrix_A[j * N + i]  * x_metod[i + j];
        }

        for (int q = 1; (q < N) && (i - q >= 0); q++) {
            nev[i] -= matrix_A[(i - q) * N + q] * x_metod[i - q];
        }

    }

    double nev_l_1 = 0.0; 
    double nev_l_2 = 0.0;
    double nev_l_inf = 0.0;

    double b_l_1 = 0.0;
    double b_l_2 = 0.0;
    double b_l_inf = 0.0;

    for (int w = 0; w < M; w++) { //

        double summ = nev[w];

        double b = vector_b[w];

        nev_l_1 += abs(summ );  // l1 
        nev_l_2 += pow(summ , 2);  // l2
        if (abs(summ ) > nev_l_inf) { // l_inf
            nev_l_inf = abs(summ );
        }


        b_l_1 += abs(b);
        b_l_2 += pow(b, 2);
        if (abs(b) > b_l_inf) { // l_inf
            b_l_inf = abs(b);
        }
    }


    nev_l_2 = sqrt(nev_l_2);
    b_l_2 = sqrt(b_l_2);




    fprintf(errors, "\n");

    fprintf(errors, "|                                    ВЫЧИСЛЕНИЕ НЕВЯЗОК ДЛЯ %s                                           |  \n", name);    
    fprintf(errors, "|                            |       ||*||_1          |          ||*||_2         |        ||*||_inf         |\n");
    fprintf(errors, "|  абсоллютная               |    %14.8e      |     %14.8e       |     %14.8e       | \n", nev_l_1, nev_l_2, nev_l_inf);
    fprintf(errors, "| относительно правой части  |    %14.8e      |     %14.8e       |     %14.8e       | \n", nev_l_1 / b_l_1, nev_l_2 / b_l_2, nev_l_inf / b_l_inf);

    delete[] nev;
}


int main()
{

    double a = -3.0;
    double b = 3.0; // выбираем границы отрезка

    const int K = 5; // выбираем количество интервалов

    const int N = 5; // 

    const int M = K * (N - 1) + 1; // общее количество узлов M интерполяции (красных точек)

    const int L = 10; // количесвто случайных точек на каждом из интервалов

    const int M_viz = 300; //  количество точек по которым будем сторить график
    int M_V = (M_viz - 1) / K; // количество среднее серых точек среди двух красных(узлов) 
    M_V *= K;
    M_V += 1;
    //cout << "M_V = " << M_V << "\n";


    double h = (b - a) / (M - 1); // шаг равномерной сетки (расстояние между красными точками)

    double* znam; // знаменатели
    znam = new double[N];

    double* XX; // равномерная сетка  по узлам
    XX = new double[M];

    double* F_XX; // равномерная сетка
    F_XX = new double[M];

    double t = a;
    for (int i = 0; i < M; i++) {
        XX[i] = t;
        F_XX[i] = my_func(t);
        t += h;
    }

    double* SS; // точrb по которым будем сторить график
    SS = new double[M_V];

    double* F_SS; // точки для отрисовки искомого графика 
    F_SS = new double[M_V];


    // точки для отрисовки
    t = a;
    double st = (b - a) / (M_V - 1);
    for (int i = 0; i < M_V; i++) { // заполним эти точки
        SS[i] = t;
        F_SS[i] = my_func(t);
        t += st;
    }

    double* y_tilda_SS; // приближение в них (по которым будем сторить график)
    y_tilda_SS = new double[M_V];



    // так как узлы на равномерной сетке, то знаменатели будут повторяться на каждом конечном элементе(подотрезке)
    // вычислим знаменатель на каждом из подотрезков


    int pr;
    for (int j = 0; j < N; j++) { // Номер начальной точки на отрезке
        pr = 1;
        for (int n = 0; n < N; n++) {
            if (n != j) {
                pr *= (j - n);
            }
        }
        znam[j] = pow(h, N - 1) * pr;
    }


    // массив случанйх точек на всём [a,b]  их L штук на каждом из K подотрезков
    double random_t[K * L];
    for (int k = 0; k < K; k++) {


        std::uniform_real_distribution <double> distr{ XX[k * (N - 1)],  XX[(k + 1) * (N - 1)] };
        std::mt19937 gen{ std::random_device().operator ()() };

        for (int l = 0; l < L; l++) {

            random_t[k * L + l] = distr(gen);



        }
    }




    // создадим необходимые матрицы

    //cout << "Step 1" << "\n";

    VectorXd Vector_b = VectorXd::Zero(M);

    MatrixXd  Matrix_Lenta = MatrixXd::Zero(M, N); // тут столбцы = ненулевые диагонали большой матрицы( от главной и выше) 
                                                   //их всего штук N - размер стороны маленького квадратика

    for (int w = 0; w < K; w++) { // проходимся по всем интервалам
        // XX[w*(N-1)]
        // XX[(w+1) * (N - 1)]

        for (int k = w * (N - 1); k <= (w + 1) * (N - 1); k++) { // проходимся по узлам внтури w-ого интервала( то есть по строкам квадратика маленького)



            for (int n = k; n <= (w + 1) * (N - 1); n++) { //  цикл по столбцам маленького квадратика ( n  начинаем считать от k(главной диагонали) в силу симметричности)

                double summ_a = 0.0;


                for (int i = w * L; i < (w + 1) * L; i++) { //  цикл по случайным точкам внутри w-Ого интервала

                    double chisl_1 = 1.0;
                    double chisl_2 = 1.0;
                    double delta;

                    for (int p = w * (N - 1); p <= (w + 1) * (N - 1); p++) { //  цикл по узловым(красным) точкам на данном подотрезке для большого числителя 

                        delta = (random_t[i] - XX[p]);
                        if (p != k) { //  числитель для phi_k(x_i)

                            chisl_1 *= delta;
                        }

                        if (p != n) { //  числитель получаем для phi_n(x_i)

                            chisl_2 *= delta;
                        }

                    }
                    double phi_k_x_i = chisl_1 / znam[k - w * (N - 1)];
                    double phi_n_x_i = chisl_2 / znam[n - w * (N - 1)];

                    summ_a += phi_k_x_i * phi_n_x_i; // суммируем чтобы получить в итоге a_kn
                }

                Matrix_Lenta(k, n - k) += summ_a; // случай когда n = k это главная диагональ, в ленточной матрице это нулевой столбец

            }

            double summ_b = 0.0; // это b_k
            for (int i = w * L; i < (w + 1) * L; i++) { //  цикл по случайным точкам внутри w-Ого интервала

                double chisl_1 = 1.0;
                double delta;

                for (int p = w * (N - 1); p <= (w + 1) * (N - 1); p++) { //  цикл по узловым(красным) точкам на данном подотрезке для большого числителя 

                    delta = (random_t[i] - XX[p]);
                    if (p != k) { //  числитель для phi_k(x_i)

                        chisl_1 *= delta;
                    }

                }
                double phi_k_x_i = chisl_1 / znam[k - w * (N - 1)];

                summ_b += my_func(random_t[i]) * phi_k_x_i;
            }

            Vector_b(k) += summ_b;
        }
    }

    //cout << "Step 2" << "\n";

    // для подачи в решатель развернём ленточную матрицу в обычную
    MatrixXd  Matrix_A = MatrixXd::Zero(M, M);

    for (int k = 0; k < M; k++) {

        for (int n = 0; n < N; n++) {

            double el = Matrix_Lenta(k, n);

            if (n + k < M) { // условие чтобы нули из ленточной не добавлялись
                Matrix_A(k, n + k) = el;

                if (n > 0) {
                    Matrix_A(n + k, k) = el;
                }
            }


        }

    }

    double BASG_b[M];

    for (int k = 0; k < M; k++) {
        BASG_b[k] = Vector_b(k);
    }

    // решение от бииблиотеки
        // РЕШАЕМ НАШУ СИСТЕМУ ЛИНЕЙНЫХ УРАВНЕНИЙ
    VectorXd Vector_C = Matrix_A.colPivHouseholderQr().solve(Vector_b); // получаем вектор  С длинной (M)

    /////////////////////////
    ///   Метод Гаусса   ////
    /////////////////////////

    double GAUSS_A[M][N]; // ленточная матрица 


    double GAUSS_x[M]; // вектор ответа
    double GAUSS_b[M]; //  вектор правой части ( будет меняться)


    //  копируем ленточную матрицу сюда  и вектор b
    for (int k = 0; k < M; k++) {

        GAUSS_b[k] = BASG_b[k];

        for (int n = 0; n < N; n++) {

            GAUSS_A[k][n] = Matrix_Lenta(k, n);
        }
    }

    // прямой ход
    for (int m = 0; m < M; m++) {

        double  el = GAUSS_A[m][0];
       
        //  i  это индекс нижележащей модифицированной строки 
        // таких строк тут должно быть N - 1
        for (int i = 1; (i + m < M) && (i < N); i++) {
            GAUSS_b[m + i] -= GAUSS_b[m] * GAUSS_A[m][i] / el;
            for (int j = 0 ; j < N - i ; j++) {
           
                 GAUSS_A[m + i][j] -= GAUSS_A[m][i + j] * GAUSS_A[m][i ] / el ;
                
            }
            
        }
    }

    //  обратный ход
    for (int k = M - 1 ; k >= 0 ; k--) {
        double sum = 0.0;

        
        for (int j = 1 ; (k + j < M ) && ( j < N) ; j++) {

            sum += GAUSS_A[k][j] * GAUSS_x[k + j ];
        }

        GAUSS_x[k] = ( GAUSS_b[k] - sum ) / GAUSS_A[k][0];
    }


    /////////////////////////
    ///   LU-разложение  ////
    /////////////////////////

    double Lo[M][N]; // она нижнетреугольная, но мы будем рассматривать её как верхнетреугольную(то есть транспонированную версию)
    double Up[M][N];
    double LU_x[M];
    double LU_y[M];



    // заполним L и U нулями и 1 где нужно
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            Up[i][j] = 0.0;
            if (j == 0) {
                Lo[i][j] = 1.0; // в матрице L на диагонали стоят 1
            }
            else {
                Lo[i][j] = 0.0;
            }
        }
    }

    // заполнение векторов нулями
    for (int i = 0; i < M; i++) {
        LU_x[i] = 0.0;
        LU_y[i] = 0.0;
    }

    for (int i = 0; i < M; i++) {

        // вычисление диагональных элеметов U[i][0]
        for (int k = 1; (k < i + 1) && (k < N) ; k++) {
                Up[i][0] -= Up[i - k][k] * Lo[i - k][k];
        }
        Up[i][0] += Matrix_Lenta(i, 0);

        for (int j = 1; (j < N); j++) {
            // вычисление внедиагоальных элементов матриц LU разложения
            for (int k = 1; (k < i + 1) && (j + k < N); k++) {
                Up[i][j] -= Up[i - k][j + k] * Lo[i - k][k];
                Lo[i][j] -= Lo[i - k][j + k] * Up[i - k][k];
            }
            Up[i][j] += Matrix_Lenta(i, j);
            Lo[i][j] += Matrix_Lenta(i, j);
            Lo[i][j] /= Up[i][0];
        }
    }


    for (int m = 0; m < M; m++) {
        double sum = 0.0;
        for (int j = 1; j < N && m - j >= 0; j++) {
            sum += LU_y[m - j] * Lo[m - j][j];
        }
        LU_y[m] = BASG_b[m] - sum;
    }

    for (int m = M - 1; m >= 0; m--) {
        double sum = 0.0;
        for (int j = 1; j < N && m + j < M; j++) {
            sum += LU_x[m + j] * Up[m][j];
        }
        LU_x[m] = (LU_y[m] - sum) / Up[m][0];
    }


    /////////////////////////
    ///  Метод Холецкого ////
    /////////////////////////

    double LL[M][N]; // она нижнетреугольная, но мы будем рассматривать её как верхнетреугольную(то есть транспонированную версию)
    double LL_x[M];
    double LL_y[M];



    //  
    for (int i = 0; i < M; i++) {
        for (int j = 0; j < N; j++) {
            LL[i][j] = 0.0;
        }
    }

    // заполнение векторов нулями
    for (int i = 0; i < M; i++) {
        LL_x[i] = 0.0;
        LL_y[i] = 0.0;
    }

    
    for (int i = 0; i < M; i++) {

        // вычисление диагональных элеметов LL[i][0]
        if (i == 0) {
            LL[i][0] = sqrt(Matrix_Lenta(i, 0));
        }
        else {
            for (int k = 1; (k < i + 1) && (k < N); k++) {
                LL[i][0] -= LL[i - k][k] * LL[i - k][k];
            }
            LL[i][0] += Matrix_Lenta(i, 0);
            LL[i][0] = sqrt(LL[i][0]);
        }

        for (int j = 1; (j < N); j++) {
            // вычисление внедиагоальных элементов матриц LL разложения
            if (i == 0) {
                LL[i][j] = Matrix_Lenta(i, j) / LL[i][0];
            }
            else {
                for (int k = 1; (k < i + 1) && (j + k < N); k++) {
                    LL[i][j] -= LL[i - k][j + k] * LL[i - k][k];
                }
                LL[i][j] += Matrix_Lenta(i, j);
                LL[i][j] /= LL[i][0];
            }
			
        }
    }


    // решаем нижнетреугольную  матрицу
    for (int m = 0; m < M; m++) {
        double sum = 0.0;
        for (int j = 1; j < N && m - j >= 0; j++) {
            sum += LL_y[m - j] * LL[m - j][j];
        }
        LL_y[m] = ( BASG_b[m] - sum ) / LL[m][0];
    }
    // решаем верхнетреугольную матрицу
    for (int m = M - 1; m >= 0; m--) {
        double sum = 0.0;
        for (int j = 1; j < N && m + j < M; j++) {
            sum += LL_x[m + j] * LL[m][j];
        }
        LL_x[m] = (LL_y[m] - sum) / LL[m][0];
    }


    ///////////////////////////////////////////////////
    // метод последовательных верхних релаксаций ПВР //
    ///////////////////////////////////////////////////

    double PVR_delta = 1e-20;
    double PVR_eps = 1e-10;
    double w = 1.5;

    double PVR_x[M];
    double PVR_xpred[M];

    // x0 начальное приближенное решение - вектор из нулей
    for (int i = 0; i < M; i++) {
        PVR_x[i] = 0.0;
        PVR_xpred[i] = 0.0;
    }

    int PVR_N = 10;
    int PVR_cur = 0;

    double PVR_res = 0.0;
    double PVR_norm_cur = 0.0;

    int iter = 0;

    while (PVR_cur < PVR_N) {
        iter++;
        // обновление вектора по координатно
        for (int i = 0; i < M; i++) {
            double tmpprev = 0.0;
            double tmpcur = 0.0;
            for (int j = 1; (j < N) && (i + j < M); j++) {
                tmpprev += Matrix_Lenta(i, j) * PVR_xpred[i + j];
            }
            for (int q = 1; (q < N) && (i - q >= 0); q++) {
                tmpcur += Matrix_Lenta(i - q, q) * PVR_x[i - q];
            }
            PVR_x[i] = (1.0 - w) * PVR_xpred[i] + w * (BASG_b[i] - tmpcur - tmpprev) / Matrix_Lenta(i, 0);
        }

        // вычисление L2 нормы z(k)-z(k-1)
        PVR_res = 0.0;
        for (int m = 0; m < M; m++) {
            PVR_res += (PVR_x[m] - PVR_xpred[m]) * (PVR_x[m] - PVR_xpred[m]);
        }
        PVR_res = sqrt(PVR_res);

        // вычисление L2 нормы z(k)
        PVR_norm_cur = 0.0;
        for (int m = 0; m < M; m++) {
            PVR_norm_cur = PVR_x[m] * PVR_x[m];
        }
        PVR_norm_cur = sqrt(PVR_norm_cur);

        // обновление предыдущего вектора приближения
        for (int m = 0; m < M; m++) {
            PVR_xpred[m] = PVR_x[m];
        }

        // если критерий удолетворен
        // увеличиваем счетчик,  то есть надо идти маленькими шажочками 10 раз подряд
        // L2 длинну шага сравниваем с долей решения
        if (PVR_res < PVR_eps * PVR_norm_cur + PVR_delta) {
            PVR_cur++;
        }
        else {
            // инчае обнуляем счетчик
            PVR_cur = 0;
        }

        // погрешность 
        double nev = 0.0;

        for (int m = 0; m < M; m++) {
            nev += (PVR_x[m] - Vector_C(m)) * (PVR_x[m] - Vector_C(m) );
        }
        nev = sqrt(nev);

       
        //cout << " PVR_cur =  " << PVR_cur << "  , PVR_res = "<< PVR_res<< " , porg =  " << nev << "\n";
    }

    cout << " PVR iterations  =  " << iter << "\n";

    ////////////////////////////////////
    // метод сопряженных градиентов СГ //
    ////////////////////////////////////

    double SG_eps = 1e-10;

    double SG_x[M], SG_r[M], SG_z[M];
    double SG_Az[M], SG_Ax[M];
    double alpha;
    double SG_res = 0.0;

    // x0 начальное приближенное решение - вектор из нулей
    for (int i = 0; i < M; i++) {
        SG_x[i] = 0.0;
    }


    // векторно-матричное умножение SG_Ax = A_ * SG_x
    for (int i = 0; i < M; i++) {
        SG_Ax[i] = 0.0;
        for (int j = 0; (j < N) && (i + j < M); j++) {
            SG_Ax[i] += Matrix_Lenta(i, j) * SG_x[i + j];
        }
        for (int q = 1; (q < N) && (i - q >= 0); q++) {
            SG_Ax[i] += Matrix_Lenta(i - q, q) * SG_x[i - q];
        }
    }

    for (int i = 0; i < M; i++) {
        // r0 = b - Ax0
        SG_r[i] = BASG_b[i] - SG_Ax[i];
        // z0 = r0
        SG_z[i] = SG_r[i];
    }


    // вычисление inf нормы невязки
    for (int i = 0; i < M; i++) {
        if (abs(SG_r[i]) > SG_res) {
            SG_res = abs(SG_r[i]);
        }
    }

    // поиск решения методом сопряженных градиентов
    iter = 0;
    while (SG_res > SG_eps) {
        iter++;
        // №1 Поиск alpha
        double r = 0.0;
        // находим (r(k-1),r(k-1))
        for (int i = 0; i < M; i++) {
            r += SG_r[i] * SG_r[i];
        }

        // матрично-векторное умножение (одно за итерацию)
        for (int i = 0; i < M; i++) {
            SG_Az[i] = 0.0;
            for (int j = 0; (j < N) && (i + j < M); j++) {
                SG_Az[i] += Matrix_Lenta(i, j) * SG_z[i + j];
            }
            for (int q = 1; (q < N) && (i - q >= 0); q++) {
                SG_Az[i] += Matrix_Lenta(i - q, q) * SG_z[i - q];
            }
        }

        // находим (Az(k-1),z(k-1))
        double az = 0.0;
        for (int i = 0; i < M; i++) {
            az += SG_Az[i] * SG_z[i];
        }
        // alpha = (r(k-1),r(k-1))/(Az(k-1),z(k-1))
        alpha = r / az;

        //№2 Уточнение решения x(k) = x(k-1) + alpha * z(k-1)
        for (int i = 0; i < M; i++) {
            SG_x[i] += alpha * SG_z[i];
        }

        //№3 Уточнение невязки r(k) = r(k-1) - alpha * A*z(k-1)
        for (int i = 0; i < M; i++) {
            SG_r[i] -= alpha * SG_Az[i];
        }

        //4. находим beta

        // вычисление (r(k),r(k))
        double rr = 0.0;
        for (int i = 0; i < M; i++) {
            rr += SG_r[i] * SG_r[i];
        }
        // beta = (r(k),r(k)) / (r(k-1),r(k-1))
        double beta = rr / r;

        //5. Уточнение z

        // z(k) = r(k) + beta * z(k-1)
        for (int i = 0; i < M; i++) {
            SG_z[i] = SG_r[i] + beta * SG_z[i];
        }


        // вычисление inf нормы невязки
        SG_res = 0.0;
        for (int i = 0; i < M; i++) {
            if (abs(SG_r[i]) > SG_res) {
                SG_res = abs(SG_r[i]);
            }
        }
    }

    cout << " SG iterations  =  " << iter << "\n";


    // ВЫЧИСЛЕНИЕ НОРМ


   FILE* errors;
    errors = fopen("errors.txt", "w");

    double vector_C[M];
    for (int w = 0; w < M; w++) { // проходимся по всем интервалам
        vector_C[w] = Vector_C(w);  }



    pogr(M, (double*) &vector_C, (double*) &GAUSS_x, "Гаусс", errors);
    pogr(M, (double*) &vector_C, (double*) &LU_x,    "LU", errors);
    pogr(M, (double*) &vector_C, (double*) &LL_x,    "Холецкий", errors);
    pogr(M, (double*) &vector_C, (double*) &PVR_x,   "ПВР", errors);
    pogr(M, (double*) &vector_C, (double*) &SG_x,    "СГ", errors);
    
    fprintf(errors, "------------------------------------------------------------------------------------------------ \n");

    double Matrixx_A[M][N]; // ленточная матрица 

    //  копируем ленточную матрицу сюда  и вектор b
    for (int k = 0; k < M; k++) {

        for (int n = 0; n < N; n++) {

            Matrixx_A[k][n] = Matrix_Lenta(k, n);
        }
    }


    nevyazka(M, N, (double*) &BASG_b, (double*) &Matrixx_A, (double*) &GAUSS_x, "Гаусс", errors);
    nevyazka(M, N, (double*) &BASG_b, (double*) &Matrixx_A, (double*) &LU_x,    "LU", errors);
    nevyazka(M, N, (double*) &BASG_b, (double*) &Matrixx_A, (double*) &LL_x,    "Холецкий", errors);
    nevyazka(M, N, (double*) &BASG_b, (double*) &Matrixx_A, (double*) &PVR_x,   "ПВР", errors);
    nevyazka(M, N, (double*) &BASG_b, (double*) &Matrixx_A, (double*) &SG_x,    "СГ", errors);

    fprintf(errors, "------------------------------------------------------------------------------------------------------------- \n");

    fclose(errors);
    

    cout << " Work is done ! ";

    return 0;

}


