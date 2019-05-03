// Source code for Greek letter risk estimation of American-Asian options.
// Line 117 returns Greek estimations via Likelihood Ratio Method (LRM).
// Line 118 returns Greek estimations via Pathwise Estimation (PE).
// As submitted, the program output returns Greek estimations via LRM.
// Comment out line 117 nad convert line 118 to code to have output return Greek estimations via PE.

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <random>

using namespace std;
using namespace std::chrono;

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);

double expiration_time = 1.0;
double initial_stock_price = 50;
double strike_price = 55;
double dividend_yield = 0;
double risk_free_rate = 0.05;
double sigma = 0.3;
double no_of_divisions = 100;

double delta_T = expiration_time / ((double)no_of_divisions);
double delta_R = (risk_free_rate - 0.5*pow(sigma, 2))*delta_T;
double delta_SD = sigma*sqrt(delta_T);

double get_normal()
{
    std::normal_distribution <double> distribution(0,1);
    double number = distribution(generator);
    return (number);
}

void Greeks_LRM(double stock_price) {
    double delta_lrm = 0;
    double vega_lrm = 0;
    double gamma_lrm = 0;
    //double mean;
    int N_trials = 1000000;
    
    for (int i =0; i<N_trials; i++)
    {
        vector<double>z(no_of_divisions,0);
        double mean = stock_price;
        double vega_summer = 0;
        vector<double>current_stock_price(no_of_divisions,0);
        //vector<double>d(no_of_divisions,0);
        //vector<double>d_diff(no_of_divisions,0);
        current_stock_price[0] = stock_price;
        for (int j = 1; j<no_of_divisions; j++)
        {
            z[j] = get_normal();
            //cout << z[j] << endl;
            current_stock_price[j] = current_stock_price[j-1]*exp(delta_R + delta_SD * z[j]);
            mean += current_stock_price[j];
            vega_summer += (((pow(z[j],2) - 1)/sigma) - (z[j]*sqrt(delta_T)));
        }
        mean = mean/no_of_divisions;
        delta_lrm += ((exp(-risk_free_rate*expiration_time))*max((strike_price - mean),0.0))*z[1]/(stock_price*delta_SD);
        gamma_lrm += ((exp(-risk_free_rate*expiration_time))*max((strike_price - mean),0.0))*(((pow(z[1], 2) - 1)/((pow(stock_price, 2))*(pow(delta_SD, 2)))) - (z[1]/((pow(stock_price, 2))*delta_SD)));
        double vega_full = vega_summer*(exp(-risk_free_rate*expiration_time))*max((strike_price - mean),0.0);
        vega_lrm += vega_full;
    }
    delta_lrm = delta_lrm/N_trials;
    gamma_lrm = gamma_lrm/N_trials;
    vega_lrm = vega_lrm/N_trials;
    cout << stock_price << "\t" << delta_lrm << "\t" << gamma_lrm << "\t" << vega_lrm << endl;
}

void Greeks_Pathwise(double stock_price) {
    double delta_pe = 0;
    double vega_pe = 0;
    //double mean;
    int N_trials = 1000000;
    
    for (int i =0; i<N_trials; i++)
    {
        vector<double>z(no_of_divisions,0);
        double mean = stock_price;
        //double vega_summer = 0;
        vector<double>current_stock_price(no_of_divisions,0);
        //vector<double>d(no_of_divisions,0);
        //vector<double>d_diff(no_of_divisions,0);
        current_stock_price[0] = stock_price;
        for (int j = 1; j<no_of_divisions; j++)
        {
            z[j] = get_normal();
            //cout << z[j] << endl;
            current_stock_price[j] = current_stock_price[j-1]*exp(delta_R + delta_SD * z[j]);
            mean += current_stock_price[j];
        }
        mean = mean/no_of_divisions;
        if(mean > strike_price)
        {
            delta_pe += exp(-risk_free_rate*expiration_time)*(mean/stock_price)  ;
            vega_pe += exp(-risk_free_rate*expiration_time) * (log(mean/stock_price) - (risk_free_rate + 0.5*pow(sigma, 2))*expiration_time) * (mean/sigma);
        }
        
    }
    delta_pe = delta_pe/N_trials;
    vega_pe = vega_pe/N_trials;
    cout << stock_price << "\t" << delta_pe << "\t" << vega_pe << endl;
}


int main(int argc, char* argv[])
{
    double S = 100;
    for (double i=1; i<=S; i++) {
        Greeks_LRM(i);
        //Greeks_Pathwise(i);
    }
    return 0;
};
