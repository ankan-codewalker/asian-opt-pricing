// Source code for Longstaff-Schwartz pricing of American-Asian options.
// Citation - Professor Ramavarapu Sreenivas' code samples.

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <random>
#include "/Users/post2ankan/Downloads/newmat10/newmat.h"

using namespace std;
using namespace std::chrono;

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine generator(seed);

double initial_stock_price = 50;
double risk_free_rate = 0.05;
double strike_price = 55;
double volatility = 0.3;
double expiration_time = 1.0;
double seconds;
int no_of_simulations = 100000;
int no_of_divisions = 1000;
int no_of_trials;
time_t t1, t2;

double delta_T = expiration_time / ((double)no_of_divisions);
double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_T;
double delta_SD = volatility*sqrt(delta_T);
double R = exp(risk_free_rate*expiration_time / ((double)no_of_divisions));


double max(double a, double b) {
    return (b < a) ? a : b;
}

int min(int a, int b) {
    return (b < a) ? b : a;
}

float get_uniform()
{
    //    return (((float) random())/(pow(2.0, 31.0)-1.0));
    uniform_real_distribution<double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return number;
}

double get_gaussian()
{
    return (sqrt(-2.0*log(get_uniform()))*cos(2 * 3.141592654*get_uniform()));
}

double get_normal()
{
    std::normal_distribution <double> distribution(0,1);
    double number = distribution(generator);
    return (number);
}


Matrix polynomial_regression(Matrix Independent_Variables, Matrix Dependent_Variable, int order, int no_of_observations)
{
    Matrix X(no_of_observations, order);
    Matrix Y(no_of_observations, 1);
    
    for (int i = 1; i <= no_of_observations; i++)
        Y(i, 1) = Dependent_Variable(i, 1);
    
    for (int j = 1; j <= order; j++)
        for (int i = 1; i <= no_of_observations; i++)
            X(i, j) = pow(Independent_Variables(i, 1), j - 1);
    
    // return inv(XT*X)*XT*Y
    Matrix X_transpose_times_X(order, order);
    X_transpose_times_X = X.t()*X;
    return (X_transpose_times_X.i() * X.t() * Y);
}

double Longstaff_Schwartz_Put(int no_of_simulations) {
    double Asian_Put_Price = 0.0;
    for (int k = 0; k < (no_of_simulations / 200 + 1); k++) {
        
        if (k != no_of_simulations / 200)
            no_of_trials = 200;
        else
            no_of_trials = no_of_simulations % 200;
        
        if (no_of_trials != 0) {
            double stock_price[200][1000];
            double mean_stock_price[200][1000];
            
            for (int i = 0; i < no_of_trials; i++)
            {
                stock_price[i][0] = initial_stock_price;
                mean_stock_price[i][0] = initial_stock_price;
            }
            
            for (int i = 0; i < no_of_trials; i++)
                for (int j = 1; j < no_of_divisions; j++)
                {
                    stock_price[i][j] = stock_price[i][j - 1] * exp(delta_R + delta_SD*get_gaussian());
                    mean_stock_price[i][j] = ((j)*mean_stock_price[i][j - 1] + stock_price[i][j]) / (j + 1);
                }
            
            
            
            double value[no_of_trials];
            // initialize the value based on the price at final stage
            for (int i = 0; i < no_of_trials; i++)
                value[i] = max(0.0, strike_price - mean_stock_price[i][no_of_divisions - 1]);
            
            for (int i = (no_of_divisions - 1); i > 0; i--) {
                Matrix independent_variables(no_of_trials, 1);
                Matrix dependent_variables(no_of_trials, 1);
                int no_of_variables = 0;
                for (int j = 0; j < no_of_trials; j++) {
                    if (max(0.0, strike_price - mean_stock_price[j][i]) > 0) {
                        no_of_variables++;
                        independent_variables(no_of_variables, 1) = mean_stock_price[j][i];   // Should it be just stock_price???
                        dependent_variables(no_of_variables, 1) = value[j] / R;
                    }
                }
                
                if (no_of_variables > 0)
                {
                    
                    // regressing the dependent_variables on the independent variables using a 4th order polynomial
                    Matrix a(min(5, no_of_variables), 1);
                    a = polynomial_regression(independent_variables, dependent_variables, min(5, no_of_variables), no_of_variables);
                    if (no_of_variables >= 5) {
                        for (int j = 0; j < no_of_trials; j++) {
                            if (((strike_price - mean_stock_price[j][i]) >(a(1, 1) + (a(2, 1)*mean_stock_price[j][i]) + (a(3, 1)*pow(mean_stock_price[j][i], 2)) + (a(4, 1)*pow(mean_stock_price[j][i], 3)) +
                                                                           (a(5, 1)*pow(mean_stock_price[j][i], 4)))) &&
                                ((strike_price - mean_stock_price[j][i]) > 0.0))
                                value[j] = strike_price - mean_stock_price[j][i];
                            else
                                value[j] = value[j] / R;
                        }
                    }
                    else if (no_of_variables == 4) {
                        for (int j = 0; j < no_of_trials; j++) {
                            if (((strike_price - mean_stock_price[j][i]) >(a(1, 1) + (a(2, 1)*mean_stock_price[j][i]) + (a(3, 1)*pow(mean_stock_price[j][i], 2)) + (a(4, 1)*pow(mean_stock_price[j][i], 3)))) &&
                                ((strike_price - mean_stock_price[j][i]) > 0.0))
                                value[j] = strike_price - mean_stock_price[j][i];
                            else
                                value[j] = value[j] / R;
                        }
                        
                    }
                    else if (no_of_variables == 3) {
                        for (int j = 0; j < no_of_trials; j++) {
                            if (((strike_price - mean_stock_price[j][i]) >(a(1, 1) + (a(2, 1)*mean_stock_price[j][i]) + (a(3, 1)*pow(mean_stock_price[j][i], 2)))) &&
                                ((strike_price - mean_stock_price[j][i]) > 0.0))
                                value[j] = strike_price - mean_stock_price[j][i];
                            else
                                value[j] = value[j] / R;
                        }
                    }
                    else if (no_of_variables == 2) {
                        for (int j = 0; j < no_of_trials; j++) {
                            if (((strike_price - mean_stock_price[j][i]) >(a(1, 1) + (a(2, 1)*mean_stock_price[j][i]))) &&
                                ((strike_price - mean_stock_price[j][i]) > 0.0))
                                value[j] = strike_price - mean_stock_price[j][i];
                            else
                                value[j] = value[j] / R;
                        }
                    }
                    else {
                        for (int j = 0; j < no_of_trials; j++) {
                            if (((strike_price - mean_stock_price[j][i]) > a(1, 1)) &&
                                ((strike_price - mean_stock_price[j][i]) > 0.0))
                                value[j] = strike_price - mean_stock_price[j][i];
                            else
                                value[j] = value[j] / R;
                        }
                    }
                    
                }
            }
            double local_put_price = 0.0;
            for (int j = 0; j < no_of_trials; j++)
                local_put_price += value[j];
            local_put_price = (local_put_price / ((float)no_of_trials)) / R;
            Asian_Put_Price += local_put_price;
        }
    }
    return (Asian_Put_Price / ((double)no_of_simulations / 200));
}

double Longstaff_Schwartz_Call(int no_of_simulations) {
    double Asian_Call_Price = 0.0;
    for (int k = 0; k < (no_of_simulations / 200 + 1); k++) {
        
        if (k != no_of_simulations / 200)
            no_of_trials = 200;
        else
            no_of_trials = no_of_simulations % 200;
        
        if (no_of_trials != 0) {
            double stock_price[200][1000];
            double mean_stock_price[200][1000];
            
            for (int i = 0; i < no_of_trials; i++)
            {
                stock_price[i][0] = initial_stock_price;
                mean_stock_price[i][0] = initial_stock_price;
            }
            
            for (int i = 0; i < no_of_trials; i++)
                for (int j = 1; j < no_of_divisions; j++)
                {
                    stock_price[i][j] = stock_price[i][j - 1] * exp(delta_R + delta_SD*get_gaussian());
                    mean_stock_price[i][j] = ((j)*mean_stock_price[i][j - 1] + stock_price[i][j]) / (j + 1);
                }
            
            
            
            double value[no_of_trials];
            // initialize the value based on the price at final stage
            for (int i = 0; i < no_of_trials; i++)
                value[i] = max(0.0, mean_stock_price[i][no_of_divisions - 1] - strike_price);
            
            for (int i = (no_of_divisions - 1); i > 0; i--) {
                Matrix independent_variables(no_of_trials, 1);
                Matrix dependent_variables(no_of_trials, 1);
                int no_of_variables = 0;
                for (int j = 0; j < no_of_trials; j++) {
                    if (max(0.0, mean_stock_price[j][i] - strike_price) > 0) {
                        no_of_variables++;
                        independent_variables(no_of_variables, 1) = mean_stock_price[j][i];   // Should it be just stock_price???
                        dependent_variables(no_of_variables, 1) = value[j] / R;
                    }
                }
                
                if (no_of_variables > 0)
                {
                    
                    // regressing the dependent_variables on the independent variables using a 4th order polynomial
                    Matrix a(min(5, no_of_variables), 1);
                    a = polynomial_regression(independent_variables, dependent_variables, min(5, no_of_variables), no_of_variables);
                    if (no_of_variables >= 5) {
                        for (int j = 0; j < no_of_trials; j++) {
                            if (((mean_stock_price[j][i] - strike_price) >(a(1, 1) + (a(2, 1)*mean_stock_price[j][i]) + (a(3, 1)*pow(mean_stock_price[j][i], 2)) + (a(4, 1)*pow(mean_stock_price[j][i], 3)) +
                                                                           (a(5, 1)*pow(mean_stock_price[j][i], 4)))) &&
                                ((mean_stock_price[j][i] - strike_price) > 0.0))
                                value[j] = mean_stock_price[j][i] - strike_price;
                            else
                                value[j] = value[j] / R;
                        }
                    }
                    else if (no_of_variables == 4) {
                        for (int j = 0; j < no_of_trials; j++) {
                            if (((mean_stock_price[j][i] - strike_price) >(a(1, 1) + (a(2, 1)*mean_stock_price[j][i]) + (a(3, 1)*pow(mean_stock_price[j][i], 2)) + (a(4, 1)*pow(mean_stock_price[j][i], 3)))) &&
                                ((mean_stock_price[j][i] - strike_price) > 0.0))
                                value[j] = mean_stock_price[j][i] - strike_price;
                            else
                                value[j] = value[j] / R;
                        }
                        
                    }
                    else if (no_of_variables == 3) {
                        for (int j = 0; j < no_of_trials; j++) {
                            if (((mean_stock_price[j][i] - strike_price) >(a(1, 1) + (a(2, 1)*mean_stock_price[j][i]) + (a(3, 1)*pow(mean_stock_price[j][i], 2)))) &&
                                ((mean_stock_price[j][i] - strike_price) > 0.0))
                                value[j] = mean_stock_price[j][i] - strike_price;
                            else
                                value[j] = value[j] / R;
                        }
                    }
                    else if (no_of_variables == 2) {
                        for (int j = 0; j < no_of_trials; j++) {
                            if (((mean_stock_price[j][i] - strike_price) >(a(1, 1) + (a(2, 1)*mean_stock_price[j][i]))) &&
                                ((mean_stock_price[j][i] - strike_price) > 0.0))
                                value[j] = mean_stock_price[j][i] - strike_price;
                            else
                                value[j] = value[j] / R;
                        }
                    }
                    else {
                        for (int j = 0; j < no_of_trials; j++) {
                            if (((mean_stock_price[j][i] - strike_price) > a(1, 1)) &&
                                ((mean_stock_price[j][i] - strike_price) > 0.0))
                                value[j] = mean_stock_price[j][i] - strike_price;
                            else
                                value[j] = value[j] / R;
                        }
                    }
                    
                }
            }
            double local_call_price = 0.0;
            for (int j = 0; j < no_of_trials; j++)
                local_call_price += value[j];
            local_call_price = (local_call_price / ((float)no_of_trials)) / R;
            Asian_Call_Price += local_call_price;
        }
    }
    return (Asian_Call_Price / ((double)no_of_simulations / 200));
}


int main(int argc, char* argv[])
{
    int N = 100000;
    cout << "--------------------------------" << endl;
    cout << "American Asian Call Option Pricing using Longstaff and Schwartz's Least Squares Monte Carlo Simulation" << endl;
    cout << "Expiration Time (Years) = " << expiration_time << endl;
    cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
    cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
    cout << "Initial Stock Price = " << initial_stock_price << endl;
    cout << "Strike Price = " << strike_price << endl;
    cout << "Number of Simulations = " << no_of_simulations << endl;
    cout << "Number of Divisions = " << no_of_divisions << endl;
    cout << "R = " << R << endl;
    cout << "--------------------------------" << endl;
    cout << "Stock" << "\t" << "Put" << "\t" << "Time_P" << "\t" << "Call" << "\t" << "Time_C" << endl;
    
    for (int i = 100; i < N; i = i+100) {
        auto start_put = high_resolution_clock::now();
        Longstaff_Schwartz_Put(i);
        auto stop_put = high_resolution_clock::now();
        auto duration_put = duration_cast<milliseconds>(stop_put - start_put);
        auto start_call = high_resolution_clock::now();
        Longstaff_Schwartz_Call(i);
        auto stop_call = high_resolution_clock::now();
        auto duration_call = duration_cast<milliseconds>(stop_call - start_call);
        cout << i << "\t" << Longstaff_Schwartz_Put(i) << "\t" << ((double)duration_put.count()/1000) << "\t" << Longstaff_Schwartz_Call(i) << "\t" << ((double)duration_call.count()/1000) << endl;
    }
}
