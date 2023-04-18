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
#include <vector>
#include <map>
#include <tuple>
#include <stdlib.h>

using namespace std;

int n=20, m=30, d=5*360, d1=2*360;
double R=exp(0.05*1.0/360), sig=0.2, u=exp(sig*sqrt(1.0/360)), p=(exp(R*1.0/360)-1/u)/(u-1/u), q = 1-p, c=4;
double triggerprice=120, strikeprice=110, conversionrat=1, S0=100;
int B = floor(log(triggerprice/S0) / log(u));

unsigned int ***M;
int ***S;
double ***h;
double ***P;
double **P2;

//map <tuple<int, int>, int> M;

unsigned int nChoosek( int k,  int n ){
    if (k > n || k < 0)
        return 1;
    if (k==0 || n==0)
        return 1;
    if (k * 2 > n)
        k = n-k;
    if (k == 0)
        return 1;

    unsigned int result = n;
    for(unsigned int i = 2; i <= k; i++ ) {
        result *= (unsigned)(n-i+1);
        result /= (unsigned)i;
    }
    return result;
}

void construct_S(int b){
    int b_ = b-(m-1);
    for (int i=0; i<m; i++){
        for (int j=0; j<=i; j++){
            S[b][i][j] = 0;
            if (j-(i-j)>b_)
                S[b][i][j]= 1;
        }
    }

}

//void update_M(int b, unsigned int**** M){
//    int** S;
//    S = construct_S(b-(m-1));
//
//    for (int m_=1; m_<m; m_++){
//        for (int l=0; l<=m_+1; l++){
//            for (int k=0; k<=m_; k++){
//                if (l==0){
//                    if (S[m_][k]==0){
//                        if (k==m_){ M[k][b][l][m_] = M[k][b][l][m_-1]; }
//                        else {M[k][b][l][m_] = M[k+1][b][l][m_-1] + M[k][b][l][m_-1]; }
//                    }
//                    continue;
//                }
//
//                if (k==m_)
//                    M[k][b][l][m_] = (S[m_][k]==0) ? (M[k][b][l][m_-1]) : (M[k][b][l-1][m_-1]);
//                else
//                    M[k][b][l][m_] = (S[m_][k]==0) ? (M[k][b][l][m_-1] + M[k+1][b][l][m_-1]) : (M[k][b][l-1][m_-1] + M[k+1][b][l-1][m_-1]);
//
//            }
//        }
//    }
//
//    for (int i=0; i<m; i++){
//        delete [] S[i];
//    }
//    delete[] S;
//}


//start from update_M(0, 0, 0, b) for each b
void update_M(int i, int j, int l, int b){

    if (S[b][i][j]==1) // which means the stock price at the current node is higher than the trigger price
        l++;
    if (i==m-1){
        M[j][b][l]++;
        return;
    }

    update_M(i+1, j, l, b);
    update_M(i+1, j+1, l, b);
    
}


double calculate_h(int b, int l, int i, int j) {
    if (i<m)
        return 1;
    if (b>=m-1)
        return 1;
    else if (b<=-m)
        return 0;
    
    int b_ = b + (m - 1);  //b_ is the index of dimention b
    if (b_<0)
        return 0;
    if (b_>=m)
        return 1;
    unsigned long long sum=0;
    for (int k = 0; k <= m - 1; k++) {
        sum += M[k][b_][l]  * nChoosek(j + k - m + 1, i - m + 1);
    }
    if (sum==0)
        return 0;
    
    unsigned long long sum1=0;
    for (int k = 0; k <= b_ / 2; k++){
        sum1 += M[k][b_][l]  * nChoosek(j + k - m + 1, i - m + 1);

    }
    if (sum1 / sum < 0 || sum1 / sum > 1) { cout << "Wrong probabilities!" << endl;}
    return 1.0*sum1 / sum;
}

//double** stockpriceArray(int d) {
//    double** S = new double* [d];
//    for (int i = 0; i < d; i++) {
//        S[i] = new double[d];
//        for (int j = 0; j <= i; j++) {
//            S[i][j] = 0;
//        }
//    }
//    S[0][0] = S0;
//    for (int i = 1; i < d; i++) {
//        for (int j = 0; j <= i; j++) {
//            S[i][j] = 100 * pow(u, 2*j-i);
//        }
//    }
//
//    return S;
//}

vector<int> coupon;

void coupondays() {
    int j=1;
    double n;
    for (int i = 0; i <= d; i++) {
        double t = (5.0 / d) * i;
        double next = (5.0 / d) * (i + 1);
        n = j * 0.5;
        if (t < n && next >= n) {
            coupon.push_back(i);
            j++;
        }
    }
}

void cov_tree_next_three_years(){
    map<int, double> Stock_price;
    double S_cur = S0;
    Stock_price[0] = S0;
    for (int j = 1; j < d; j++){
        S_cur = S_cur*u;
        Stock_price[j] = S_cur;
    }
    S_cur = S0;
    for (int j=-1; j > -d; j--){
        S_cur = S_cur/u;
        Stock_price[j] = S_cur;
    }
    
    for (int i=d-1; i>=d1; i--){
        for (int j=0; j<=i; j++){
            if (i==d-1)
                P2[i-d1][j] = max(strikeprice, Stock_price[2*j-i]);
            else{
                if (Stock_price[2*j-i] > strikeprice)
                    P2[i-d1][j] = max(strikeprice, Stock_price[2*j-i] * conversionrat);
                else{
                    double contval = p * P2[i+1-d1][j+1] + q * P2[i+1-d1][j];
                    if (find(coupon.begin(), coupon.end(), i) != coupon.end())
                        contval += c;
                
                    P2[i-d1][j] = max(contval / R, Stock_price[2*j-i] * conversionrat);
                }
            }
        }
    }
    
}

double conv_tree(int i,int j, int l, double stock_price) {
    
    
    
    //cout << i << j << l << int(stock_price) << endl;
    double contval = 0;
    double pi1 = 0, pi2 = 0;
    
    if (i < d1) {
        if (i < m) {
            if (l >= n) {
                P[i][j][l] = max(strikeprice, stock_price * conversionrat);
                return P[i][j][l];
            }
            
            if (stock_price * u > triggerprice) {
                if (P[i + 1][j + 1][l + 1] < 0)
                    P[i + 1][j + 1][l + 1] = conv_tree(i + 1, j + 1, l + 1, stock_price * u);
                contval += p * P[i + 1][j + 1][l + 1];
            }
            else {
                if (P[i + 1][j + 1][l] < 0)
                    P[i + 1][j + 1][l] = conv_tree(i + 1, j + 1, l, stock_price * u);
                contval += p * P[i + 1][j + 1][l];
            }
            if (stock_price / u > triggerprice) {
                if (P[i + 1][j][l + 1] < 0)
                    P[i + 1][j][l + 1] = conv_tree(i + 1, j, l + 1, stock_price / u);
                contval += q * P[i][j][l + 1];
            }
            else {
                if (P[i + 1][j][l] < 0)
                    P[i + 1][j][l] = conv_tree(i + 1, j, l, stock_price / u);
                contval += q * P[i + 1][j][l];
            }
//double bond_price = max(0.0, stock_price - strikeprice) + strikeprice/(1+R*(d/360-i/360));
            P[i][j][l] = max(contval / R, stock_price * conversionrat);
            return P[i][j][l];
            

        }
        else {
            double prob;

            if (stock_price > triggerprice * pow(u, m)) {
                P[i][j][l] = max(strikeprice, stock_price * conversionrat);
                return P[i][j][l];
            }


            else if (stock_price < triggerprice * pow(1 / u, m)) {
                if (P[i + 1][j + 1][l] < 0)
                    P[i + 1][j + 1][l] = conv_tree(i + 1, j + 1, l, stock_price * u);
                if (P[i + 1][j][l] < 0)
                    P[i + 1][j][l] = conv_tree(i + 1, j, l, stock_price / u);
                contval = p * P[i + 1][j + 1][l] + q * P[i + 1][j][l];
                //double bond_price = max(0.0, stock_price - strikeprice) + strikeprice/(1+R*(d/360-i/360));
                P[i][j][l] = max(contval / R, stock_price * conversionrat);
                return P[i][j][l];
            }

            else {
                if (l >= n) {
                    P[i][j][l] = max(strikeprice, stock_price * conversionrat);
                    return P[i][j][l];
                }
                
                
                if (h[i][j][l]<0)
                    h[i][j][l] = calculate_h(B-(2*j-i),l,i,j);
                prob = h[i][j][l];
                if (stock_price * u > triggerprice) {
                    if (l == m) {
                        if (P[i + 1][j + 1][l] < 0)
                            P[i + 1][j + 1][l] = conv_tree(i + 1, j + 1, l, stock_price * u);
                        pi1 = P[i + 1][j + 1][l];
                    }
                    else {
                        if (P[i + 1][j + 1][l + 1] < 0)
                            P[i + 1][j + 1][l + 1] = conv_tree(i + 1, j + 1, l + 1, stock_price * u);
                        if (P[i + 1][j + 1][l] < 0)
                            P[i + 1][j + 1][l] = conv_tree(i + 1, j + 1, l, stock_price * u);
                        pi1 = prob * P[i + 1][j + 1][l + 1] + (1 - prob) * P[i + 1][j + 1][l];
                    }
                }
                else {
                    if (l == 0) {
                        if (P[i + 1][j + 1][l] < 0)
                            P[i + 1][j + 1][l] = conv_tree(i + 1, j + 1, l, stock_price * u);
                        pi1 = P[i + 1][j + 1][l];

                    }
                    else {
                        if (P[i + 1][j + 1][l] < 0)
                            P[i + 1][j + 1][l] = conv_tree(i + 1, j + 1, l, stock_price * u);
                        if (P[i + 1][j + 1][l - 1] < 0)
                            P[i + 1][j + 1][l - 1] = conv_tree(i + 1, j + 1, l - 1, stock_price * u);
                        pi1 = prob * P[i + 1][j + 1][l] + (1 - prob) * P[i + 1][j + 1][l - 1];
                    }
                }
                if (stock_price / u > triggerprice) {
                    if (l == m) {
                        if (P[i + 1][j][l] < 0)
                            P[i + 1][j][l] = conv_tree(i + 1, j, l, stock_price / u);
                        pi2 = P[i + 1][j][l];
                    }
                    else {
                        if (P[i + 1][j][l + 1] < 0)
                            P[i + 1][j][l + 1] = conv_tree(i + 1, j, l + 1, stock_price / u);
                        if (P[i + 1][j][l] < 0)
                            P[i + 1][j][l] = conv_tree(i + 1, j, l, stock_price / u);
                        pi2 += prob * P[i + 1][j][l + 1] + (1 - prob) * P[i + 1][j][l];
                    }
                }
                else {
                    if (l == 0) {
                        if (P[i + 1][j][l] < 0)
                            P[i + 1][j][l] = conv_tree(i + 1, j, l, stock_price / u);
                        pi2 = P[i + 1][j][l];
                    }
                    else {
                        if (P[i + 1][j][l] < 0)
                            P[i + 1][j][l] = conv_tree(i + 1, j, l, stock_price / u);
                        if (P[i + 1][j][l - 1] < 0)
                            P[i + 1][j][l - 1] = conv_tree(i + 1, j, l - 1, stock_price / u);
                        pi2 += prob * P[i + 1][j][l] + (1 - prob) * P[i + 1][j][l - 1];
                    }
                }
                contval = p * pi1 + q * pi2;
                if (find(coupon.begin(), coupon.end(), i) != coupon.end())
                    contval += c;
                
                    //double bond_price = max(0.0, stock_price - strikeprice) + strikeprice/(1+R*(d/360-i/360));
                P[i][j][l] = max(contval / R, stock_price * conversionrat);
                return P[i][j][l];
                
            }
        }
 
    }
    
    else{ // (i >= 360*2 )
        return P2[i-d1][j];
    }
    
    return 0;
}


int main(int argc, const char * argv[]) {
    p=0.504;
    q=1-p;
    
    //cout << nChoosek(3, 5) << " " << nChoosek(7, 5) << " " << nChoosek(200, 500) << endl;
    
//    double ** S;
//    S = new double* [d];
//    for (int i=0; i<d; i++){
//        S[i] = new double [d];
//        for (int j=0; j<d; j++)
//            S[i][j] = 0;
//    }
//
//    for (int i=0; i<d; i++)
//        for (int j=0; j<=i; j++)
//            S[i][j] = S0 * pow(u, j) * pow(1.0/u, i-j);
//
//    for (int i=0; i<d; i++){
//        for (int j=0; j<d; j++)
//            cout << S[i][j] << " ";
//        cout << endl;
//    }
//
//    cout << endl;
    
    
//    unsigned int ****M_; //k, b, l, m
//    M_ = new unsigned int*** [m];  // k from 0 to m-1
//    for (int k=0;  k < m; k++){
//        M_[k] = new unsigned int** [2*m-2];// b from -(m-1) to m-2
//        for (int b=0; b < 2*m-2; b++){
//            M_[k][b] = new unsigned int* [m+1]; // l from 0 to m
//            for (int l=0; l<=m; l++){
//                M_[k][b][l] = new unsigned int [m]; // m_ from 0 to m-1
//                for (int m_=0; m_<m; m_++){
//                    M_[k][b][l][m_] = 0;
//                }
//            }
//        }
//
//    }
//    M_[0][0+m-1][0][0] = 1;
//    M_[0][-1+m-1][1][0] = 1;
    
    coupondays();

    S = new int**[2*m-2]; // b from -(m-1) to m-2
    for (int b=0; b<2*m-2; b++){
        S[b] = new int*[m]; // i from 0 to m
        for (int i=0; i<m; i++){
            S[b][i] = new int [m]; // j from 0 to m
        }
    }
    for (int b=0; b<2*m-2; b++){
        construct_S(b);
    }

    M = new unsigned int** [m];  // k from 0 to m-1
    for (int k=0;  k < m; k++){
        M[k] = new unsigned int* [2*m-2];// b from -(m-1) to m-2
        for (int b=0; b < 2*m-2; b++){
            M[k][b] = new unsigned int [m+1]; // l from 0 to m
            for (int l=0; l<=m; l++){
                M[k][b][l] = 0;
            
            }
        }
        
    }
    for (int b=0; b<=2*m-3; b++){
        update_M(0,0,0, b);
    }
    
//    for (int k=0;  k < m; k++){
//        for (int b=0; b < 2*m-2; b++){
//            for (int l=0; l<=m; l++){
//                if (M[k][b][l]>0)
//                    cout << "M[" << k << "][" << b-(m-1) << "][" << l << "][" << m << "] = " << M[k][b][l] << endl;
//            }
//        }
//    }

    
//    for (int k=0;  k < m; k++){
//        for (int b=0; b < 2*m-2; b++){
//            for (int l=0; l<=m; l++){
//                delete[] M_[k][b][l];
//            }
//            delete[] M_[k][b];
//        }
//        delete[] M_[k];
//    }
//    delete[] M_;
    
    
    
    h = new double** [d1];  // i from 0 to d-1
    for (int i=0;  i < d1; i++){
        h[i] = new double* [d1];// j from 0 to i
        for (int j=0; j < d1; j++){
            h[i][j] = new double [m+1]; // l from 0 to m
            for (int l=0; l<=m; l++){
                h[i][j][l] = -1;

            }
        }

    }
    
//    int B, b_cur;
//
//    //cout << (pow(u, B)*S0 <= triggerprice && pow(u, B+1)*S0 > triggerprice) << endl;
//
//    for (int i=0;  i < d; i++){
//        for (int j=0; j <= i; j++){
//            b_cur = B-(2*j-i);
//            for (int l=0; l<=m; l++){
//                h[i][j][l] = calculate_h(b_cur, l ,M, i, j);
//            }
//        }
//    }
    
    P = new double** [d1+1]; // i from 0 to 2*360-1 +1
    for (int i=0; i<d1+1; i++){
        P[i] = new double* [d1+1]; // j from 0 to 2*360-1 +1
        for (int j=0; j<d1+1; j++){
            P[i][j] = new double [m+1];
            for (int l=0; l<=m; l++){ // l from 0 to m
                P[i][j][l] = -1.0;
            }
        }
    }
    
    P2 = new double* [d-d1]; // i from 2*360 to 5*360-1
    for (int i=0; i<d-d1; i++)
        P2[i] = new double [d]; //j from 0 to 5*360-1
    
    cov_tree_next_three_years();
    

    
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
    price = conv_tree(0,0,0, S0);
    cout << price << endl;
    cout << endl;
    
//    cout << "P:" << endl;
    
//    for (int l=0; l<=m; l++){
//        cout << "l=" << l << endl;
//        for (int i=0; i<d; i++){
//            for (int j=0; j<d; j++)
//                cout << P[i][j][l] << " ";
//            cout << endl;
//        }
//        cout << endl;
//    }
    
//    cout << "h:" << endl;
//    for (int l=0; l<=m; l++){
//        cout << "l=" << l << endl;
//        for (int i=0; i<d; i++){
//            for (int j=0; j<d; j++)
//                cout << h[i][j][l] << " ";
//            cout << endl;
//        }
//        cout << endl;
//    }


    for (int k=0;  k < m; k++){
        for (int b=0; b < 2*m-2; b++)
            delete[] M[k][b] ;
        delete[] M[k];
    }
    delete[] M;

    for (int i=0;  i < d1; i++){
        for (int j=0; j < d1; j++)
            delete[] h[i][j];
        delete[] h[i];
    }
    delete[] h;
    
    for (int i=0; i<d1+1; i++){
        for (int j=0; j<d1+1; j++)
            delete[] P[i][j];
        delete[] P[i];
    }
    delete[] P;
    
    for (int i=0; i<d-d1; i++)
        delete[] P2[i];
    delete[] P2;
    
    for (int b=0; b<2*m-2; b++){
        for (int i=0; i<m; i++)
            delete[] S[b][i];
        delete[] S[b];
    }
    delete[] S;
    
    return 0;
}




