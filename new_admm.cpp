#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>
#include <bitset>

using namespace std;
using namespace arma;
using namespace Rcpp;


#define ARMA_DONT_PRINT_ERRORS


void InitializeEta1(mat Z, vec Y, mat IsQ, mat VtV, mat &Eta1)
{
    const int max_iter = 30;
    const int num_var = Z.n_cols;
    const int num_sample = Z.n_rows;
    mat Eta1_old = Eta1;
    int iter = 1;
	  double aa = 0.0;
    mat I = eye<mat>(num_sample, num_sample);    
    while (iter < max_iter)
    {
        vec Y_tmp = zeros<vec>(num_sample);
        for (size_t i = 0; i < num_var; i++)
        {
            Y_tmp += Z.col(i) % Eta1.col(i);
        }

        for (size_t i = 0; i < num_var; i++)
        {
            vec Y_iter = Y - ( Y_tmp - Z.col(i) % Eta1.col(i) );
            mat Zdiag = diagmat(Z.col(i));
            Eta1.col(i) = inv(Zdiag.t() * IsQ * Zdiag + 0.1 * I) * (Zdiag * IsQ * Y_iter);
        }
		
        // to check stop criteria
        if (norm(Eta1_old - Eta1, 2) < 1e-2)
        {
            aa = norm(Eta1_old - Eta1, 2);
			      cout << "del=" << aa << endl;
			      break;
        }
        else
        {
            Eta1_old = Eta1;
        }
        iter++;
    } // end while
    cout << "iter" << iter << endl;
} // end funcs

// Date: 2019-1-7 17:49:53
// [[Rcpp::export]]
SEXP subG_ADMM_new(SEXP Yin, SEXP Xin, SEXP Zin, SEXP Eta1in, SEXP Utin, SEXP lambdain, SEXP admm_iterin, SEXP tolin, SEXP ridgein, SEXP kappain)
{
    try
    {                         // *admm Algorithm
        vec Y = as<vec>(Yin); // *a vector response, dim = num_sample x 1
        vec X = as<vec>(Xin); // *a vector we are interested, dim = num_sample x 1
        mat Z = as<mat>(Zin); // *dim num_sample x q1
        // mat W = as<mat>(Win); // covariate matrix, num_sample x c

        int admm_iter = Rcpp::as<int>(admm_iterin);
        double lambda = Rcpp::as<double>(lambdain);
        double tol = Rcpp::as<double>(tolin);
        size_t num_sample = Y.n_elem;
        size_t num_var = Z.n_cols;
        int ridge = Rcpp::as<int>(ridgein);
        double kappa = Rcpp::as<double>(kappain);

        //================================================
        // compute the Q using Z0 in the manuscript
        // Q -- num_sample x num_sample
        // mat Q = W * inv_sympd(W.t() * W) * W.t();
        
        //================================================
        // compute Z -- num_sample x num_var
        Z = Z.each_col() % X;
        //mat I = eye<mat>(num_sample, num_sample);
        // *Z_d -- num_sample x q2*num_sample
        //mat Z_d = RowWiseKronProd(I, Z);

        //=============================================
        // compute V -- num_sample*(num_sample - 1)/2 x num_sample
        mat I = eye<mat>(num_sample, num_sample);
        mat V = zeros<mat>(num_sample * (num_sample - 1) * 0.5, num_sample);
        size_t t = 0;
        for (size_t i = 0; i < (num_sample - 1); i++)
        {
            for (size_t j = i + 1; j < num_sample; j++)
            {
                V.row(t) = I.row(i) - I.row(j);
                t++;
            } // end for
        }     // end for
        // num_sample x num_sample
        mat VtV = V.t() * V;
		
        //mat I2 = eye<mat>(num_var, num_var);
        //mat v_star = KronProd(I2, V);
        //======================================================
        // Compute Eta1, initialize Eta1 -- num_sample x num_var

        //mat Eta1 = zeros<mat>(num_sample, num_var);
        mat Eta1 = as<mat>(Eta1in);
        mat IsQ = I;
        mat Ut = as<mat>(Utin);
        if (ridge >= 1)
        {
            InitializeEta1(Z, Y, IsQ, VtV, Eta1);
        }
        //========================================
        // Compute alpha
        vec Y_tmp = zeros<vec>(num_sample);
        // for (size_t i = 0; i < num_var; i++)
        // {
        //     Y_tmp += Z.col(i) % Eta1.col(i);
        // } // end for
        // vec alpha = inv(W.t() * W) * W.t() * (Y - Y_tmp);

        //================================================
        // main loop
        vec rho = ones<vec>(2) / 1; // LU CHANGED, originally, vec rho = ones<vec>(2);
        mat Theta = zeros<mat>(num_var, num_sample * (num_sample - 1) * 0.5);
        mat T = zeros<mat>(num_var, num_sample * (num_sample - 1) * 0.5);
        mat T_old = T;
        size_t iter = 0;
        double indicator = 0.0;
        double gamma = 3.0; // default in most papers
        while (iter < admm_iter)
        {
            // cout << "admm_step: iter = " << iter + 1 << endl;
            //=======================================================
            // *update theta
            // B -- num_var x num_sample*(num_sample-1)/2
            //cout << "admm_step: update theta " << endl;
            mat B = Eta1.t() * V.t() - T / rho(0);
            for (size_t i = 0; i < B.n_cols; i++)
            {
                double sign = 1.0;
                if (norm(B.col(i), 2) > lambda * gamma)
                {
                    Theta.col(i) = B.col(i);
                }
                else
                {
                    indicator = 1 - lambda * (1.0 / rho(0)) / norm(B.col(i), 2);
                    if (indicator > 0) // LU CHANGED, originally, indicator > 0
                    {
                        Theta.col(i) = indicator * Theta.col(i)/(1 - 1.0/gamma/rho(0) );
                    }
                    else
                    {
                        Theta.col(i).zeros();
                    } // end fi
                }     // end fi
            }         // end func

            //cout<<"Theta = "<<Theta.t()<<endl;

            //=====================================
            // *update Eta1
            //cout << "admm_step: update Eta1 " << endl;
            Y_tmp.zeros();
            for (size_t i = 0; i < num_var; i++)
            {
                Y_tmp += Z.col(i) % Eta1.col(i);
            } // end for

            for (size_t i = 0; i < num_var; i++)
            {
                vec Y_iter = Y - ( Y_tmp - Z.col(i) % Eta1.col(i) );
                mat Zdiag = diagmat(Z.col(i));
                vec Tcol = conv_to<vec>::from(T.row(i));
                vec Thetacol = conv_to<vec>::from(Theta.row(i));
                Eta1.col(i) = inv(Zdiag.t() * IsQ * Zdiag + rho(0) * VtV + kappa * I) * (Zdiag * IsQ * Y_iter + V.t() * Tcol + rho(0) * V.t() * Thetacol + kappa * Ut.col(i));
            } // end for
            // beta = inv(Z_d.t() * (eye(size(Q)) - Q) * Z_d + rho(0) * VtV) * (Z_d.t() * (eye(size(Q)) - Q) * Y + v_tide.t() * T + rho(0) * v_tide.t() * Theta);
            //cout<<"Sigma="<<Sigma.t()<<endl;

            //===========================================
            // *update alpha
            //cout << "admm_step: update alpha " << endl;
            // Y_tmp.zeros();
            // for (size_t i = 0; i < num_var; i++)
            // {
            //     Y_tmp += Z.col(i) % Eta1.col(i);
            // } // end for
            // alpha = inv(W.t() * W) * W.t() * (Y - Y_tmp);
            //cout<<"alpha="<<alpha.t()<<endl;
            //===========================================
            // *update T
            //cout << "admm_step: update T " << endl;
            T = T_old + rho(0) * (Theta - Eta1.t() * V.t());
            //cout<<"T="<<T.t()<<endl;

            // *stop rule
            // hist_gradT[iter] = norm(gradT, 2);
            // cout << "T = " << norm(T - T_old, 2) << endl;
            // if(norm(T - T_old, 2) < 1e-3){break;}
            if(norm(T - T_old, 2) < tol){break;}
            T_old = T;
            iter++;
            //if(norm(gradT, 2) < tol){break;}
        } // *end for admm_iter
        // cout << "iter" << iter << endl;
        return List::create(Named("Eta1") = Eta1, Named("T") = T, Named("Theta") = Theta, Named("num_iter") = iter - 1);
    }
    catch (std::exception &ex)
    {
        forward_exception_to_r(ex);
    }
    catch (...)
    {
        ::Rf_error("C++ exception (unknown reason)...");
    }
    return R_NilValue;
} // end func


/////////////////////////////////////////////////////////////////////////////////////////
//                             CODE END HERE                                           //
/////////////////////////////////////////////////////////////////////////////////////////
