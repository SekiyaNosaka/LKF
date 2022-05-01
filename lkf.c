/*
 * @author: nosaka
 * 確率ロボティクスのカルマンフィルタ実装
 * なるべく教科書内の文字や語録を参照して書いている
*/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

double randNormal(double mu, double sigma)
{
    double uniform = ((double)rand()+1.0)/((double)RAND_MAX+2.0);
    double z = sqrt(-2.0*log(uniform)) * sin(2.0*M_PI*uniform);
    return mu + sigma*z;
}


int main()
{
    srand(time(NULL));

    int dt = 1; // time step
    int measure_time = 0;
    int end_time = 50;
    
    float R = 1.0; // convariance R (where is a scalar's variance)
    float Q = 2.0; // convariance Q (where is a scalar's variance)

    float epsilon = 0.0;
    float delta = 0.0;
    
    float A = 1.2; // A in the equation of state (where is a scalar's variance)
    float B = 0.0; // B in the equation of state (where is a scalar's variance)
    float C = 1.0; // C in the measure equation (where is a scalar's variance)
    
    float x_true = 0.0;
    float x_t = 0.0;
    float x_old = 0.0;

    float l_sig_t = 0.0;
    float l_sig_old = 0.0;

    float z = 0.0;

    float K = 0.0;

    float mes_diff = 0.0;
    float mes_add = 0.0;
    float cor_diff = 0.0;
    float cor_add = 0.0;

    while (measure_time < end_time)
    {
        // 真値(実際ロボットで搭載する際は不明)
        epsilon = randNormal(0.0, R*R);
        x_true = (A*x_true) + epsilon;
        // -> 真値append point

        // 観測
        delta = randNormal(0.0, Q*Q);
        z = C*x_true + delta;
        // -> 観測append point
        mes_diff = (x_true - z);

        // 予測
        epsilon = randNormal(0.0, R*R);
        x_t = (A*x_old) + epsilon; // 予測値
        l_sig_t = A*l_sig_old*A + R; // 予測誤差分散

        // 最適カルマンゲイン計算
        K = (l_sig_t*C) / (C*l_sig_t*C + Q);

        // 計測更新
        x_old = x_t + K*(z - C*x_t); // 更新処理
        l_sig_old = (1 - K*C)*l_sig_t; // 更新された誤差の分散
        // -> 更新append point
        cor_diff = (x_true - x_old);

        // タイムカウント
        measure_time += dt;

        mes_add += (mes_diff*mes_diff);
        cor_add += (cor_diff*cor_diff);
    }

    printf("観測二乗和誤差:  %lf\n", mes_add/50);
    printf("更新二乗和誤差:  %lf\n", cor_add/50);

    return 0;
}
