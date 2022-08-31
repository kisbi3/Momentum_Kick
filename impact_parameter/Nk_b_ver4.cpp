#include <fstream>
#include <string>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <cstring>
#include <thread>
#define PI 3.1415926535897932385
#define mbp 0.938272046 // proton mass
#define mb 0.13957018 // pion mass

#include <iomanip>

//Range
float beta0 = 0.;
float betan = PI/2.;
float l0 = 0.; // jet trajectory
float ln = 0.8;
float b00 = 0.; // jet source point
float b0n = 0.8;
float phis0 = 0.;
float phisn = PI/2.;
float bi = 0.;
float bf = 1.58;
int n = 100;
//Constants
float RA = 0.8;
float RB = 0.8;
float sigma = 0.14; //0.14 11
int K = 21; //21 367
float t0 = 0.6;
float zeta = 0.2;
//1
float TATB_p(float b0, float phis, float b, float l, float beta)
{
    float squareA = pow(b0,2)+pow(l,2)+pow(b,2)/4.-b*l*cos(phis)-b*b0*cos(beta)
                            +2.*l*b0*(cos(beta)*cos(phis)+sin(beta)*sin(phis));
    float squareB = pow(b0,2)+pow(l,2)+pow(b,2)/4.+b*l*cos(phis)+b*b0*cos(beta)
                            +2.*l*b0*(cos(beta)*cos(phis)+sin(beta)*sin(phis));
    float TA_p = (3./(2.*PI*pow(RA,3)))*sqrt(pow(RA,2)-squareA);
    float TB_p = (3./(2.*PI*pow(RB,3)))*sqrt(pow(RB,2)-squareB);
    if(RA - sqrt(squareA) > 0. && RB - sqrt(squareB) > 0.)
    {return(TA_p + TB_p);}
    else return(0.);
}
float beta_Nk(float b0, float phis, float b, float l)
{
    float dbeta = (betan-beta0)/n;
    float beta, sum = 0.;
    for(int i = 0; i < n; i++)
    {
        beta = beta0 + i*dbeta;
        sum += TATB_p(b0,phis,b,l,beta)*dbeta;
    }
    return((sigma*K/(2.*(t0+l)))*sum);
}
//2
float Nk(float b0, float phis, float b)
{
    float dl = (ln-l0)/1000.;
    float l, sum = 0.;
    for(int i = 0; i < 1000; i++)
    {
        l = l0 + i*dl;
        sum += beta_Nk(b0,phis,b,l)*dl;
    }
    return(sum);
}
//3
float beta_TATB_0(float b0, float b, float beta)
{
    float squareA = pow(b0,2)+pow(b,2)/4.-b*b0*cos(beta);
    float squareB = pow(b0,2)+pow(b,2)/4.+b*b0*cos(beta);
    float TA_0 = (3./(2.*PI*pow(RA,3)))*sqrt(pow(RA,2)-squareA);
    float TB_0 = (3./(2.*PI*pow(RB,3)))*sqrt(pow(RB,2)-squareB);
    if(RA - sqrt(squareA) > 0. && RB - sqrt(squareB) > 0.)
    {return(TA_0*TB_0);}
    else return(0.);
}
float TATB_0(float b0, float b)
{
    float dbeta = (betan-beta0)/n;
    float beta, sum = 0.;
    for(int i = 0; i < n; i++)
    {
        beta = beta0 + i*dbeta;
        sum += beta_TATB_0(b0,b,beta)*dbeta;
    }
    return(sum);
}
//4
float Ib0_Deno_TATB(float b)
{
    float db0 = (b0n-b00)/n;
    float b0, sum = 0.;
    for(int i = 0; i < n; i++)
    {
        b0 = b00 + i*db0;
        sum += b0*TATB_0(b0,b)*db0;
    }
    return(sum);
}
//5
float Probability(float b0, float b)
{
    float result;
    if(Ib0_Deno_TATB(b) > 0.0001)
    {result = TATB_0(b0,b)/Ib0_Deno_TATB(b);}
    else result = 0.;
    return(result);
}
//6
float Aver_Nk(float phis, float b)
{
    float db0 = (b0n-b00)/n;
    float b0, sum1, sum2 = 0.;
    // std::cout<<std::endl<<"phis : "<<std::setw(15)<<phis<<std::endl;
    for(int i = 0; i < n; i++)
    {
        b0 = b00 + i*db0;
        // std::cout<<"b0 : "<<std::setw(15)<<b0<<std::setw(15)<<"Nk : "<<Nk(b0,phis,b)<<std::endl;
        sum1 += b0*Nk(b0,phis,b)*exp(-zeta*Nk(b0,phis,b))*Probability(b0,b)*db0;
        sum2 += b0*exp(-zeta*Nk(b0,phis,b))*Probability(b0,b)*db0;
    }
    return(sum1/sum2);
}
//7
float Iphis_Nk(float b)
{
    float dphis = (phisn-phis0)/n;
    float phis, sum = 0.;
    float dist = 0.;
    // std::cout<<std::endl<<"impact parameter : "<<b<<std::endl;
    for(int i = 0; i < n; i++)
    {
        phis = phis0 + i*dphis;
        dist = Aver_Nk(phis,b);
        sum += dist*dphis;
        std::cout<<"phis : "<<phis<<std::setw(15)<<"Nk : "<<dist<<std::setw(15)<<"sum : "<<sum<<std::endl;
    }
    float Iphis_Nk = sum;
    return((2./PI)*Iphis_Nk);
}
void Nk1()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b1.csv",std::ios::out);
    for(float b = 0.; b <= 0.09; b += 0.1)
    {
        float Nk = Iphis_Nk(b);
        fs << b << "," << Nk << std::endl;
    }
    // float b = 0.;
    // std::cout << b << "," << Nk << std::endl;
    // fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk2()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b2.csv",std::ios::out);
    // for(float b = 0.1; b < 0.2; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 0.1;
    std::cout << b << "," << Nk << std::endl;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk3()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b3.csv",std::ios::out);
    // for(float b = 0.2; b < 0.3; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 0.2;
    std::cout << b << "," << Nk << std::endl;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk4()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b4.csv",std::ios::out);
    // for(float b = 0.3; b < 0.4; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 0.3;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk5()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b5.csv",std::ios::out);
    // for(float b = 0.4; b < 0.5; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 0.4;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk6()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b6.csv",std::ios::out);
    // for(float b = 0.5; b < 0.6; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 0.5;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk7()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b7.csv",std::ios::out);
    // for(float b = 0.6; b < 0.7; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 0.6;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk8()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b8.csv",std::ios::out);
    // for(float b = 0.7; b < 0.8; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 0.7;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk9()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b9.csv",std::ios::out);
    // for(float b = 0.8; b < 0.9; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 0.8;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk10()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b10.csv",std::ios::out);
    // for(float b = 0.9; b < 1.; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 0.9;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk11()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b11.csv",std::ios::out);
    // for(float b = 1.; b < 1.1; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 1.;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk12()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b12.csv",std::ios::out);
    // for(float b = 1.1; b < 1.2; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 1.1;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk13()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b13.csv",std::ios::out);
    // for(float b = 1.2; b < 1.3; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 1.2;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk14()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b14.csv",std::ios::out);
    // for(float b = 1.3; b < 1.4; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 1.3;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk15()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b15.csv",std::ios::out);
    // for(float b = 1.4; b < 1.5; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 1.4;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
void Nk16()
{
    std::string str_buf;
    std::fstream fs;
    fs.open("Nk_b16.csv",std::ios::out);
    // for(float b = 1.5; b < 1.58; b += 0.1)
    // {
    //     float Nk = Iphis_Nk(b);
    //     fs << b << "," << Nk << std::endl;
    // }
    float b = 1.5;
    fs << b << "," << Nk << std::endl;
    fs.close();
}
int main()
{
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
    // std::thread t15(Nk15);
    // std::thread t16(Nk16);
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
    // t15.join();
    // t16.join();

    // std::string filename = "cat Nk_b1.csv Nk_b2.csv Nk_b3.csv Nk_b4.csv Nk_b5.csv Nk_b6.csv Nk_b7.csv Nk_b8.csv Nk_b9.csv Nk_b10.csv Nk_b11.csv Nk_b12.csv Nk_b13.csv Nk_b14.csv Nk_b15.csv Nk_b16.csv > Nk.csv";
    std::string filename = "cat Nk_b1.csv Nk_b2.csv Nk_b3.csv > Nk.csv";
    system(filename.c_str());

    return 0;
}