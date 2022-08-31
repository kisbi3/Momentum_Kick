// 분리안함
#include <fstream>
#include <string>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <thread>
#define PI 3.1415926535897932385
//===================================================================================================
void Nk1();
// void Nk2();
// void Nk3();
// void Nk4();
// void Nk5();
// void Nk6();
// void Nk7();
// void Nk8();
// void Nk9();
// void Nk10();
// void Nk11();
// void Nk12();
// void Nk13();
// void Nk14();
//===================================================================================================
double R = 0.8; //R, b0, L
double aM = PI/2.;
int n = 10;
int N = 10;

// double sigma = 0.14; // jet-medium-parton scattering cross section
double sigma = 11.;
double kappa = 367.; // Nmp/Nparti
double t0 = 0.43;
double cons = (3.*sigma*kappa)/(4.*PI*pow(R,3));
double zeta = 0.2; // jet attenuation coefficient

double dc = R  / n;
double da = aM / n;
//===================================================================================================
double sumPJET(double b)
{
    double b0x = 0.;
    double b0y = 0.;
    double dv  = R / N;
    double Z   = 0.;
    double G   = 0.;
    double sum = 0.;

    for(int i = 1; i <= N; i++)
    {
        b0x = i * dv;
        for(int j = 1; j <= N; j++)
        {
            b0y = j * dv;
            Z = R*R - b0x*b0x - b0y*b0y - b*b/4.;
            G = b0x*b;
            if(Z > G)
            {sum += sqrt(Z*Z - G*G) * dv * dv;}
            else break;
        }
    }
    return(sum);
}
double AVER2NK(double b)
{
    double ps  = 0.;
    double b0x = 0.;
    double b0y = 0.;
    
    double l = 0.;
    double C = 0.;
    double D = 0.;
    double Z = 0.;
    double G = 0.;

    double sum1 = 0.;
    double sum2 = 0.;
    double nk   = 0.;
    double EXP  = 0.;
    double pjet = 0.;

    double sumNK   = 0.;
    double AVER1NK = 0.;

    double D_PJET = sumPJET(b);

    // for(int i = 1; i <= n; i++)
    //          integral phi 0 to pi/2
    for(int i = n; i >= 1; i--)
    {
        ps = i * da;
        // ps = 0.2;
        sum1 = 0.;
        sum2 = 0.;

        for(int j = 1; j <= n; j++)
        {
            b0x = j * dc;
            for(int k = 1; k <= n; k++)
            {
                b0y = k * dc;

// Line 96 ~ 101 : integrate b0
    
                // std::cout << ps << "," << b0x << "," << b0y << std::endl;
                sumNK = 0.;
                for(int s = 1; s <= n; s++)
                {
                    l = s * dc;
                    C = R*R - (b0x*b0x + b0y*b0y + l*l + b*b/4. + 2.*l*(b0x*cos(ps) + b0y*sin(ps)));
                    D = b * (l*cos(ps) + b0x);
                    if(C > D)
                    {sumNK += (sqrt(C+D) + sqrt(C-D))/(t0+l) * dc;}
                    else break;
                }
                nk = cons * sumNK;

                // std::cout << ps << "," << b0x << "," << b0y << "," << nk << std::endl;
                if(nk != 0)
                {
                    EXP = exp(-zeta * nk);
                    std::cout<<b<<"\t"<<nk<<"\t"<<EXP<<std::endl;
                    Z = R*R - b0x*b0x - b0y*b0y - b*b/4.;
                    G = b0x*b;
                    if(Z > G)
                    {
                        pjet = sqrt(Z*Z - G*G) / D_PJET;
                        // std::cout << b0x << "," << b0y << "," << pjet << std::endl;
                        // AVER1NK numerator & denominator
                        sum1 += nk * EXP * pjet * dc * dc;
                        sum2 +=      EXP * pjet * dc * dc;
                    }
                    else break;
                }
                else break;
                // std::cout << b0x << "," << b0y << "," << pjet << std::endl;
            }
        }
        AVER1NK += sum1/sum2 * da;
        // std::cout << ps << "," << sum1 << "," << sum2 << "," << AVER1NK << std::endl;
    }
    return(AVER1NK/aM);
}
int main()
{
    std::chrono::steady_clock sc;   // create an object of `steady_clock` class
    auto start = sc.now();          // start timer
    // std::cout << AVER2NK(1.2) << std::endl;
    std::thread t1(Nk1);
    // std::thread t2(Nk2);
    // std::thread t3(Nk3);
    // std::thread t4(Nk4);
    // std::thread t5(Nk5);
    // std::thread t6(Nk6);
    // std::thread t7(Nk7);
    // std::thread t8(Nk8);
    // std::thread t9(Nk9);
    // std::thread t10(Nk10);
    // std::thread t11(Nk11);
    // std::thread t12(Nk12);
    // std::thread t13(Nk13);
    // std::thread t14(Nk14);
    t1.join();
    // t2.join();
    // t3.join();
    // t4.join();
    // t5.join();
    // t6.join();
    // t7.join();
    // t8.join();
    // t9.join();
    // t10.join();
    // t11.join();
    // t12.join();
    // t13.join();
    // t14.join();

    // std::string filename1 = "cat Nk_b1.csv Nk_b2.csv Nk_b3.csv Nk_b4.csv Nk_b5.csv Nk_b6.csv Nk_b7.csv Nk_b8.csv Nk_b9.csv Nk_b10.csv Nk_b11.csv Nk_b12.csv Nk_b13.csv Nk_b14.csv > Nk_c_res.csv";
    // system(filename1.c_str());
    // std::string filename2 = "rm Nk_b1.csv Nk_b2.csv Nk_b3.csv Nk_b4.csv Nk_b5.csv Nk_b6.csv Nk_b7.csv Nk_b8.csv Nk_b9.csv Nk_b10.csv Nk_b11.csv Nk_b12.csv Nk_b13.csv Nk_b14.csv";
    // system(filename2.c_str());

    std::string filename1 = "cat Nk_b1.csv Nk_b2.csv > Nk_c_res.csv";
    system(filename1.c_str());
    // std::string filename2 = "rm Nk_b1.csv Nk_b2.csv Nk_b3.csv Nk_b4.csv Nk_b5.csv Nk_b6.csv Nk_b7.csv Nk_b8.csv Nk_b9.csv Nk_b10.csv Nk_b11.csv";
    // system(filename2.c_str());

    auto end = sc.now();       // end timer (starting & ending is done by measuring the time at the moment the process started & ended respectively)
    auto time_span = static_cast<std::chrono::duration<double>>(end - start);   // measure time span b0ytween start & end
    std::cout << "Operation took: " << time_span.count() << " seconds" << std::endl;

    return 0;
}
// void Nk1()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b1.csv", std::ios::out);
//     // fs << 0 << "," << AVER2NK(0.) << std::endl;
//     for(double b = 0.; b <= 0.11; b += 0.01)
//     {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
void Nk1()
{
    std::string str_buf;
    std::fstream fs;
    AVER2NK(0.00);
    AVER2NK(0.01);
    AVER2NK(0.02);
    // fs.open("Nk_b1.csv", std::ios::out);
    // fs << 0.54 << "," << AVER2NK(0.54) << std::endl;
    // // fs << 0.05 << "," << AVER2NK(0.05) << std::endl;
    // // for(double b = 0.04; b <= 0.11; b += 0.01)
    // // {fs << b << "," << AVER2NK(b) << std::endl;}
    // fs.close();
}
// void Nk2()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b2.csv", std::ios::out);
//     fs << 0.55 << "," << AVER2NK(0.55) << std::endl;
//     // fs << 0.07 << "," << AVER2NK(0.07) << std::endl;
//     // for(double b = 0.12; b < 0.23; b += 0.01)
//     // {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk3()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b3.csv", std::ios::out);
//     fs << 0.08 << "," << AVER2NK(0.08) << std::endl;
//     fs << 0.09 << "," << AVER2NK(0.09) << std::endl;
//     // for(double b = 0.23; b < 0.34; b += 0.01)
//     // {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk4()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b4.csv", std::ios::out);
//     fs << 0.1 << "," << AVER2NK(0.1) << std::endl;
//     fs << 0.11 << "," << AVER2NK(0.11) << std::endl;
//     // for(double b = 0.34; b < 0.45; b += 0.01)
//     // {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk5()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b5.csv", std::ios::out);
//     fs << 0.17 << "," << AVER2NK(0.17) << std::endl;
//     fs << 0.18 << "," << AVER2NK(0.18) << std::endl;
//     // for(double b = 0.45; b < 0.56; b += 0.01)
//     // {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk6()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b6.csv", std::ios::out);
//     fs << 0.19 << "," << AVER2NK(0.19) << std::endl;
//     fs << 0.2 << "," << AVER2NK(0.2) << std::endl;
//     // for(double b = 0.56; b < 0.67; b += 0.01)
//     // {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk7()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b7.csv", std::ios::out);
//     fs << 0.21 << "," << AVER2NK(0.21) << std::endl;
//     fs << 0.22 << "," << AVER2NK(0.22) << std::endl;
//     // for(double b = 0.67; b < 0.78; b += 0.01)
//     // {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk8()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b8.csv", std::ios::out);
//     fs << 0.29 << "," << AVER2NK(0.29) << std::endl;
//     fs << 0.3 << "," << AVER2NK(0.3) << std::endl;
//     fs << 0.31 << "," << AVER2NK(0.31) << std::endl;
//     // for(double b = 0.78; b < 0.89; b += 0.01)
//     // {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk9()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b9.csv", std::ios::out);
//     fs << 0.32 << "," << AVER2NK(0.32) << std::endl;
//     fs << 0.33 << "," << AVER2NK(0.33) << std::endl;
//     // for(double b = 0.89; b < 1.; b += 0.01)
//     // {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk10()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b10.csv", std::ios::out);
//     fs << 0.41 << "," << AVER2NK(0.41) << std::endl;
//     fs << 0.42 << "," << AVER2NK(0.42) << std::endl;
//     // for(double b = 1.; b < 1.11; b += 0.01)
//     // {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk11()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b11.csv", std::ios::out);
//     fs << 0.43 << "," << AVER2NK(0.43) << std::endl;
//     fs << 0.44 << "," << AVER2NK(0.44) << std::endl;
//     // for(double b = 1.11; b < 1.22; b += 0.01)
//     // {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk12()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b12.csv", std::ios::out);
//     // fs << 1.1 << "," << AVER2NK(1.1) << std::endl;
//     for(double b = 1.22; b < 1.33; b += 0.01)
//     {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk13()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b13.csv", std::ios::out);
//     // fs << 1.2 << "," << AVER2NK(1.2) << std::endl;
//     for(double b = 1.33; b < 1.44; b += 0.01)
//     {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }
// void Nk14()
// {
//     std::string str_buf;
//     std::fstream fs;
//     fs.open("Nk_b14.csv", std::ios::out);
//     // fs << 1.3 << "," << AVER2NK(1.3) << std::endl;
//     for(double b = 1.44; b < 1.6; b += 0.01)
//     {fs << b << "," << AVER2NK(b) << std::endl;}
//     fs.close();
// }