//
//  main.cpp
//  week6
//
//  Created by 仲崽 on 10/7/22.
//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <map>
#include <tuple>
#include <stdlib.h>

using namespace std;

int n=20, m=30, d=90;
double r=0.9917, sig=0.3, u=exp(sig*sqrt(1.0/d)), p=(exp(r*1.0/d-1/u)/(u-1/u)), q = 1-p;
double triggerprice=110, strikeprice=104, conversionrat=1;

int ***h;
int ***M;
double ***P;

//map <tuple<int, int>, int> M;

int nChoosek( int n, int k ){
    if (k > n || k < 0)
        return 0;
    if (k==0 || n==0)
        return 1;
    if (k * 2 > n)
        k = n-k;
    if (k == 0)
        return 1;

    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}

int** construct_S(int b){
    int ** S = new int* [m];
    for (int i=0; i<m; i++){
        S[i] = new int[m];
        for (int j=0; j<=i; j++){
            S[i][j] = 0;
            if (j-(i-j)>b)
                S[i][j]= 1;
        }
    }
    return S;
}

void update_M(int b, int**** M){
    int** S;
    S = construct_S(b-(m-1));
    
    for (int m_=1; m_<m; m_++){
        for (int l=0; l<=m_+1; l++){
            for (int k=0; k<=m_; k++){
                if (l==0){
                    if (k==0){ if (S[m_][k]==0){ M[k][b][l][m_] = M[k+1][b][l][m_-1]; } }
                    if (k==m-1){ if (S[m_][k]==0){ M[k][b][l][m_] = M[k-1][b][l][m_-1]; } }
                    continue;
                }
                
                if (k==0)
                    M[k][b][l][m_] = (S[m_][k]==0) ? (M[k+1][b][l][m_-1]) : (M[k+1][b][l-1][m_-1]);
                else if (k==m-1)
                    M[k][b][l][m_] = (S[m_][k]==0) ? (M[k-1][b][l][m_-1]) : (M[k-1][b][l-1][m_-1]);
                else
                    M[k][b][l][m_] = (S[m_][k]==0) ? (M[k-1][b][l][m_-1] + M[k+1][b][l][m_-1]) : (M[k-1][b][l-1][m_-1] + M[k+1][b][l-1][m_-1]);
                
            }
        }
    }
    
    for (int i=0; i<m; i++){
        delete [] S[i];
    }
    delete[] S;
}


double calculate_h(int b, int l,int*** M, int i, int j) {
    
    if (i<30)
        return 1;
    
    if (b>=29)
        return 1;
    else if (b<=-30)
        return 0;
    int b_ = b + (m - 1);
    int sum=0;
    for (int k = 0; k <= m - 1; k++) {
        sum += M[k][b_][l] * nChoosek(j + k - m + 1, i - m + 1);
    }
    if (sum > 0) {
        int sum1=0;
        for (int k = 0; k <= b; k++){
            sum1 += M[k][b_][l] * nChoosek(j + k - m + 1, i - m + 1);
        }
        return (double)(sum1) / sum;
    }
    else {
        return 0;
    }
}

double** stockpriceArray(int d) {
    double** S = new double* [d];
    for (int i = 0; i < d; i++) {
        S[i] = new double[d];
        for (int j = 0; j <= i; j++) {
            S[i][j] = 0;
        }
    }
    S[0][0] = 100;
    double u = 1.007481;
    for (int i = 1; i < d; i++) {
        for (int j = 0; j <= i; j++) {
            S[i][j] = 100 * pow(u, 2*j-i);
        }
    }

    return S;
}

double conv_tree(int i,int j, int l, double stock_price) {
    double contval = 0;
    double pi1 = 0, pi2 = 0;
    if (i == d-1) {
        return max(104.00, stock_price * conversionrat);
    }
    if (l>=m || l>i){
        return max(strikeprice, stock_price*conversionrat); // not sure here
    }
    if (i < 30) {
        if (stock_price * u >triggerprice) {
            if (P[i+1][j+1][l+1]<0)
                P[i+1][j+1][l+1] = conv_tree(i + 1, j + 1, l + 1, stock_price*u);
            contval += p * P[i+1][j+1][l+1];
        }
        else {
            if (P[i+1][j+1][l]<0)
                P[i+1][j+1][l] = conv_tree(i + 1, j + 1, l, stock_price*u);
            contval += p * P[i+1][j+1][l];
        }
        if (stock_price / u > triggerprice) {
            if (P[i][j][l+1]<0)
                P[i][j][l+1] = conv_tree(i, j, l + 1, stock_price/u);
            contval += q * P[i][j][l+1];
        }
        else {
            if (P[i][j][l]<0)
                P[i][j][l] = conv_tree(i, j, l, stock_price/u);
            contval += q * P[i][j][l];
        }

        if (l >= 20) {
            return max(strikeprice, stock_price*conversionrat);
        }
        else {
            return max(contval/r, stock_price * conversionrat);
        }
 
    }
    else {
        double dummy = stock_price;
        int b = 0;
        double prob;
        while (dummy <= triggerprice && b!=30) {
            dummy = stock_price * pow(u, b) * pow((1-u), (30 - 1 - b));
            b += 1;
        }
        if (stock_price > triggerprice * pow(u, 29))
            return max(strikeprice, stock_price * conversionrat);

        else if (stock_price < triggerprice * pow(1/u, 29)){
            if (P[i+1][j+1][l]<0)
                P[i+1][j+1][l] = conv_tree(i + 1, j + 1, l, stock_price * u);
            if (P[i+1][j][l]<0)
                P[i+1][j][l] = conv_tree(i + 1, j, l, stock_price / u);
            contval = p * P[i+1][j+1][l] + q * P[i+1][j][l];
            return max(contval / r, stock_price * conversionrat);
        }
        else {
            prob = h[i][j][l];
            if (stock_price * u > triggerprice ) {
                if (P[i+1][j+1][l+1]<0)
                    P[i+1][j+1][l+1] = conv_tree(i + 1, j + 1, l + 1, stock_price * u);
                if (P[i+1][j+1][l]<0)
                    P[i+1][j+1][l] = conv_tree(i + 1, j + 1, l, stock_price * u);
                pi1 = prob * P[i+1][j+1][l+1] + (1 - prob) * P[i+1][j+1][l];
            }
            else {
                if (P[i+1][j+1][l]<0)
                    P[i+1][j+1][l] = conv_tree(i + 1, j + 1, l, stock_price * u);
                if (P[i+1][j+1][l-1]<0)
                    P[i+1][j+1][l-1] = conv_tree(i + 1, j + 1, l - 1, stock_price * u);
                pi1 = prob * P[i+1][j+1][l] + (1 - prob) * P[i+1][j+1][l-1];
            }
            if (stock_price / u > triggerprice) {
                if (P[i+1][j][l+1]<0)
                    P[i+1][j][l+1] = conv_tree(i + 1, j, l + 1, stock_price / u);
                if (P[i+1][j][l]<0)
                    P[i+1][j][l] = conv_tree(i + 1, j, l, stock_price / u);
                pi2 += prob * P[i+1][j][l+1] + (1 - prob) * P[i+1][j][l];
            }
            else {
                if (P[i+1][j+1][l]<0)
                    P[i+1][j+1][l] = conv_tree(i + 1, j + 1, l, stock_price / u);
                if (P[i+1][j+1][l-1]<0)
                    P[i+1][j+1][l-1] = conv_tree(i + 1, j + 1, l - 1, stock_price / u);
                pi2 += prob * P[i+1][j+1][l] + (1 - prob) * P[i+1][j+1][l-1];
            }
            contval = p * pi1 + q * pi2;
            if (l >= 20) {
                return max(strikeprice, stock_price * conversionrat);
            }
            else {
                return max(contval / r, stock_price * conversionrat);
            }
        }
    }
    return 0;
}


int main(int argc, const char * argv[]) {
    
    int ****M_; //k, b, l, m
    M_ = new int*** [m];  // k from 0 to m-1
    for (int k=0;  k < m; k++){
        M_[k] = new int** [2*m-3];// b from -(m-1) to m-2
        for (int b=0; b <= 2*m-3; b++){
            M_[k][b] = new int* [m+1]; // l from 0 to m
            for (int l=0; l<=m; l++){
                M_[k][b][l] = new int [m]; // m_ from 0 to m-1
                for (int m_=0; m_<m; m_++){
                    M_[k][b][l][m_] = 0;
                }
            }
        }
        
    }
    M_[0][0+m-1][0][0] = 1;
    M_[0][-1+m-1][1][0] = 1;
    
    for (int b=0; b<=2*m-3; b++){
        update_M(b, M_);
    }
    

    M = new int** [m];  // k from 0 to m-1
    for (int k=0;  k < m; k++){
        M[k] = new int* [2*m-3];// b from -(m-1) to m-2
        for (int b=0; b <= 2*m-3; b++){
            M[k][b] = new int [m+1]; // l from 0 to m
            for (int l=0; l<=m; l++){
                M[k][b][l] = M_[k][b][l][m-1];
                
            }
        }
        
    }
    
    for (int k=0;  k < m; k++){
        for (int b=0; b < 2*m-3; b++){
            for (int l=0; l<=m; l++){
                delete[] M_[k][b][l];
            }
            delete[] M_[k][b];
        }
        delete[] M_[k];
    }
    delete[] M_;
    
    
    
    h = new int** [d];  // i from 0 to d-1
    for (int i=0;  i < d; i++){
        h[i] = new int* [d];// j from 0 to i
        for (int j=0; j < d; j++){
            h[i][j] = new int [m+1]; // l from 0 to m
            for (int l=0; l<=m; l++){
                h[i][j][l] = 0;
                
            }
        }
        
    }
    
    int B=0, b_cur;
    
    for (int i=0;  i < d; i++){
        for (int j=0; j <= i; j++){
            b_cur = B-(2*j-i);
            for (int l=0; l<=m; l++)
                h[i][j][l] = calculate_h(b_cur, l ,M, i, j);
        }
    }
    
    P = new double** [d];
    for (int i=0; i<d; i++){
        P[i] = new double* [d];
        for (int j=0; j<d; j++){
            P[i][j] = new double [m+1];
            for (int l=0; l<=m; l++){
                P[i][j][l] = -1.0;
            }
        }
    }
    
//    S = stockpriceArray(d); do not need the stock price tree
    // we suppose that all past price is lower than H
//    double S_past[m];
//    double s_cur = 100;
//    for (int i=m-1; i>=0; i--){
//        r = rand() % 2; // random number 0 or 1
//        if (r==0)
//            s_cur = s_cur / u;
//        else
//            s_cur = s_cur * u;
//        S_past[i] = s_cur;
//
//    }
    
    double price;
    price = conv_tree(0,0,0, 100);
    cout << price << endl;

    
    return 0;
}

