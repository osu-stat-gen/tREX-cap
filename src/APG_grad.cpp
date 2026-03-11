#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <cmath>
#include <float.h>
#include <math.h>
#include <omp.h>
#include <R.h>  
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

#define _USE_MATH_DEFINES
#define RCPP_ARMADILLO_RETURN_ANYVEC_AS_VECTOR
#define EPSILON DBL_EPSILON
using namespace Rcpp;
using namespace arma;
using uint = unsigned int;

// https://github.com/USCbiostats/r-parallel-benchmark
//OMP parallel + SIMD
// https://chryswoods.com/vector_c++/features.html  SIMD for dummies


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]

// [[Rcpp::export]]
arma::vec rcpp_gen_angle(){
    double x1, x2, x3;
    x1 = R::runif(0.0,2*M_PI);
    // x2 = (R::runif(0,1) <0.5)? R::runif(0.0,M_PI): R::runif(M_PI, 3*M_PI/2);
    x2 = R::runif(0.0,2*M_PI);
    // if (R::runif(0,1) < 0.5)
    //     x2 = R::runif(0.0,M_PI/2);
    // else
    //     x2 = R::runif(M_PI/2, 3*M_PI/2);
    // x2 = ifelse(runif(1,0,1)<0.5, runif(1,0.0,M_PI/2), runif(1,M_PI/2, 3*M_PI/2));
    x3 = R::runif(0.0,M_PI);
    return {x1,x2,x3};
}

///////////////////////////////////
///////Rx
///////////////////////////////////

// [[Rcpp::export]]
arma::mat Rx(double theta){
    arma::mat R(3,3);
    R(0,0) = 1.0;
    R(1,1) = cos(theta);
    R(1,2) = -sin(theta);
    R(2,1) = sin(theta);
    R(2,2) = cos(theta);
    return R;
}

///////////////////////////////////
///////Ry
///////////////////////////////////

// [[Rcpp::export]]
arma::mat Ry(double theta){
    arma::mat R(3,3);
    R(0,0) = cos(theta);
    R(0, 2) = sin(theta);
    R(1, 1) = 1.0;
    R(2, 0) = -sin(theta);
    R(2, 2) = cos(theta);
    return R;
}


///////////////////////////////////
///////Rz
///////////////////////////////////
// // [[Rcpp::export]]
// arma::mat Rz(double theta){
//     arma::mat R(3,3);
//     R(0, 0) = cos(theta);
//     R(0, 1) = -sin(theta);
//     R(1, 0) = sin(theta);
//     R(1, 1) = cos(theta);
//     R(2, 2) = 1.0;
//     return R;
// }

// [[Rcpp::export]]
arma::mat Rz(double theta){
    arma::mat R(3,3);
    R(0, 0) = cos(theta);
    R(0, 1) = sin(theta);
    R(1, 0) = -sin(theta);
    R(1, 1) = cos(theta);
    R(2, 2) = 1.0;
    return R;
}



// //rowsum, supposedly fastest
// // https://stackoverflow.com/questions/51774646/efficiency-of-matrix-rowsums-vs-colsums-in-r-vs-rcpp-vs-armadillo
// // [[Rcpp::export]]
// arma::vec Arma_rowSums(const arma::mat& x) {
//   return arma::sum(x, 1);
// }
///////////////////////////////////
///////matrix distance
///////////////////////////////////
//similar speed
// [[Rcpp::export]]
arma::mat dist_trans(const arma::vec& phi, const arma::mat& x, const arma::mat& y){
    arma::mat x_sum = arma::sum(pow(x, 2), 1);
    arma::mat y_sum = arma::sum(pow(y, 2), 1);
    // arma::mat rot1 = Rx(phi(0));
    arma::mat result = -2 * x * Rz(phi(2)) * Ry(phi(1)) * Rx(phi(0)) * trans(y);
    //arma::mat result = -2 * x * Rz(phi(2))* Ry(phi(1))* Rx(phi(0))* trans(y);
   // arma::mat result = -2*x * trans(Rx(phi(2)) * Rx(phi(1)) * Rx(phi(0)) *y);
    result.each_col() += x_sum;
    result.each_row() += reshape(y_sum, 1, y.n_rows);
    return result;
}


// new code to try out to deal with the error handling: 08/28/2025
inline bool bad(double x){ return !std::isfinite(x); }

struct RootResult {
  double root;
  double err;
  int iters;
  int status; // 0=ok, 1=no bracket, 2=nan/inf, 3=maxit, 4=bad inputs
};

inline RootResult pack(double r, double e, int it, int st){
  return {r,e,it,st};
}

// tiny “zero” and safe epsilon
static inline double macheps() { return std::numeric_limits<double>::epsilon(); }
static inline bool near_zero(double v, double tol){ return std::abs(v) <= tol; }



///////////////////////////////////
///////partial derivative for beta
///////////////////////////////////
// [[Rcpp::export]]
double beta_func(double beta, const arma::mat& llambdax, const arma::mat& y12, const arma::mat& ddd) {
    arma::mat log_d = log(ddd);
    double y = arma::accu(log_d % pow(ddd, beta) % llambdax - y12 % log_d);
    return y;
}


// this is the old code that Meng had wrote
///////////////////////////////////
///////re-writing R uniroot code for CPP, specifically for this
///////////////////////////////////
//source code:
//https://svn.r-project.org/R/trunk/src/library/stats/src/zeroin.c
arma::vec r_zeroin2(double a, double b, const arma::mat& llambdax, const arma::mat& y12, const arma::mat& ddd,
                    double tol=0.000001, int Maxit=1000) {
  double c=a;
  double fa=beta_func(a, llambdax, y12, ddd);
  double fb=beta_func(b, llambdax, y12, ddd);
  double fc=fa;
  int maxit=Maxit+1;

  /* First test if we have found a root at an endpoint */
    if(fa == 0.0) {
        return {a, 0, 0};
    }
    if(fb ==  0.0) {
        return {b, 0, 0};
    }

    if (fa*fb>0.0) {
        // need better error handling here
        //does not return where this happened when called inside NLL.
        stop("f() values at end points not of opposite sign");
        return {0,0,10};   //same sign, try again
    }

    while(maxit--)	{	/* Main iteration loop	*/
        double prev_step = b-a;
        if( fabs(fc) < fabs(fb) ){				/* Swap data for b to be the	*/
            a = b;  b = c;  c = a;	/* best approximation		*/
            fa=fb;  fb=fc;  fc=fa;
        }
        double tol_act = 2*EPSILON*fabs(b) + tol/2;
        double new_step = (c-b)/2;
        //Rprintf("Rep %i, new_step=%f, fb=%f, func=%f\n", maxit, new_step, fb, beta_func(b, llambdax,y12,ddd));
        if( fabs(new_step) <= tol_act || fb == 0.0 )
        {
            return {b, fabs(c-b), Maxit-maxit};			/* Acceptable approx. is found	*/
        }

        /* Decide if the interpolation can be tried	*/
        if( fabs(prev_step) >= tol_act && fabs(fa) > fabs(fb) ) {
            double t1,t2;
            double p,q;
            double cb = c-b;
            if( a==c ) {
                t1 = fb/fa;
                p = cb*t1;
                q = 1.0 - t1;
            }
            else {
                q = fa/fc;
                t1 = fb/fc;
                t2 = fb/fa;
                p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
                q = (q-1.0) * (t1-1.0) * (t2-1.0);
            }
            if( p>0.0 )		/* p was calculated with the */
                q = -q;			/* opposite sign; make p positive */
            else			/* and assign possible minus to	*/
                p = -p;			/* q				*/
            if( (p < (0.75*cb*q-fabs(tol_act*q)/2))  && (p < fabs(prev_step*q/2)) )
                new_step = p/q;
        }

        if( fabs(new_step) < tol_act) {
            if( new_step > 0.0 )
                new_step = tol_act;
            else
                new_step = -tol_act;
        }
        a = b;	fa = fb;			/* Save the previous approx. */
        b += new_step;	fb = beta_func(b, llambdax, y12, ddd);
        if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
            c = a;  fc = fa;
        }
    }
    /* failed! */
    warning("Algorithm did not converge");
    return {b, fabs(c-b), -1,};  //did not converge
}

// this augments the brent method with a grid search+newton if the endpoints are not of opposite sign
// [[Rcpp::export]]
// arma::vec r_zeroin2(double a, double b,
//                     const arma::mat& llambdax,
//                     const arma::mat& y12,
//                     const arma::mat& ddd,
//                     double tol = 1e-6,
//                     int Maxit = 1000,
//                     int grid_points = 50) {
//     double c = a;
//     double fa = beta_func(a, llambdax, y12, ddd);
//     double fb = beta_func(b, llambdax, y12, ddd);
//     double fc = fa;
//     int maxit = Maxit + 1;

//     // Check endpoints
//     if (fa == 0.0) {
//         double beta = (a >= 0.0 ? -1e-12 : a);
//         Rcpp::Rcout << "[r_zeroin2] Endpoint solution: beta = " << beta << "\n";
//         // arma::vec out(3); out(0)=beta; out(1)=0.0; out(2)=0.0;
//         return {beta, 0.0, 0.0};
//     }
//     if (fb == 0.0) {
//         double beta = -1e-12; // forbid 0
//         Rcpp::Rcout << "[r_zeroin2] Endpoint near-zero solution: beta = " << beta << "\n";
//         // arma::vec out(3); out(0)=beta; out(1)=0.0; out(2)=0.0;
//         return {beta, 0.0, 0.0};
//     }

//     // Fallback if no sign change
//     if (fa * fb > 0.0) {
//         Rcpp::Rcout << "[r_zeroin2] No sign change, falling back to grid+Newton...\n";

//         auto g = [&](double beta) {
//             return beta_func(beta, llambdax, y12, ddd);
//         };
//         auto gprime = [&](double beta) {
//             arma::mat log_d = log(ddd);
//             arma::mat term = pow(ddd, beta) % llambdax % log_d;
//             return arma::accu(term % log_d);
//         };

//         // Step 1: coarse grid search
//         double step = (b - a) / (grid_points - 1);
//         double beta_min = a, g_min = std::abs(g(a));
//         for (int i = 1; i < grid_points; i++) {
//             double beta = a + i * step;
//             double val = std::abs(g(beta));
//             if (val < g_min) {
//                 g_min = val;
//                 beta_min = beta;
//             }
//         }

//         // Step 2: Newton refinement
//         double beta = beta_min;
//         for (int it = 0; it < Maxit; it++) {
//             double f = g(beta);
//             double fp = gprime(beta);
//             double update = f / fp;
//             beta -= update;
            
//             if (!std::isfinite(beta)) {
//                 Rcpp::warning("Beta became NaN/Inf during Newton, aborting");
//                 // arma::vec out(3); out(0)=-1.0; out(1)=0.0; out(2)=-1.0;
//                 return {-1.0, 0.0, -1.0};
//             }

//             if (beta <= a) beta = a + 1e-12;
//             if (beta >= b) beta = -1e-12; // forbid 0

//             if (std::abs(update) < tol) {
//                 Rcpp::Rcout << "[r_zeroin2] Grid+Newton solution: beta = "
//                             << beta << " (iterations: " << it+1 << ")\n";
//                 // arma::vec out(3);
//                 // out(0) = beta;
//                 // out(1) = 0.0;
//                 // out(2) = it+1;
//                 return {beta, 0.0, it+1};
//             }
//         }

//         // Did not converge
//         Rcpp::warning("Grid+Newton did not converge");
//         if (beta >= b) beta = -1e-12; // forbid 0
//         Rcpp::Rcout << "[r_zeroin2] Grid+Newton returning last iterate: beta = "
//                     << beta << "\n";
//         // arma::vec out(3); out(0)=beta; out(1)=0.0; out(2)=-1.0;
//         // return out;
//         return {beta, 0.0, -1.0};
//     }

//     // Otherwise proceed with Brent
//     while (maxit--) {
//         double prev_step = b - a;
//         if (std::fabs(fc) < std::fabs(fb)) {
//             a = b; b = c; c = a;
//             fa = fb; fb = fc; fc = fa;
//         }
//         double tol_act = 2 * EPSILON * std::fabs(b) + tol / 2;
//         double new_step = (c - b) / 2;

//         if (std::fabs(new_step) <= tol_act || fb == 0.0) {
//             double beta = b;
//             if (beta >= 0.0) beta = -1e-12; // forbid 0
//             Rcpp::Rcout << "[r_zeroin2] Brent solution: beta = " << beta
//                         << " (iterations: " << Maxit - maxit << ")\n";
//             // arma::vec out(3); out(0)=beta; out(1)=std::fabs(c - b); out(2)=Maxit - maxit;
//             return {beta, std::fabs(c - b), Maxit - maxit};
//         }

//         // Interpolation
//         if (std::fabs(prev_step) >= tol_act && std::fabs(fa) > std::fabs(fb)) {
//             double t1, t2, p, q;
//             double cb = c - b;
//             if (a == c) {
//                 t1 = fb / fa;
//                 p = cb * t1;
//                 q = 1.0 - t1;
//             } else {
//                 q = fa / fc;
//                 t1 = fb / fc;
//                 t2 = fb / fa;
//                 p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
//                 q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
//             }
//             if (p > 0.0) q = -q; else p = -p;
//             if (p < (0.75 * cb * q - std::fabs(tol_act * q) / 2) &&
//                 p < std::fabs(prev_step * q / 2))
//                 new_step = p / q;
//         }

//         if (std::fabs(new_step) < tol_act)
//             new_step = (new_step > 0.0 ? tol_act : -tol_act);

//         a = b; fa = fb;
//         b += new_step; fb = beta_func(b, llambdax, y12, ddd);
//         if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
//             c = a; fc = fa;
//         }
//     }

//     // Brent did not converge
//     double beta = b;
//     if (beta >= 0.0) beta = -1e-12; // forbid 0
//     Rcpp::warning("Brent did not converge");
//     Rcpp::Rcout << "[r_zeroin2] Brent returning last iterate: beta = " << beta << "\n";
//     // arma::vec out(3); out(0)=beta; out(1)=std::fabs(c - b); out(2)=-1.0;
//     return {beta, std::fabs(c - b), -1.0};
// }


// ver 1
// arma::vec r_zeroin2(double a, double b,
//                     const arma::mat& llambdax,
//                     const arma::mat& y12,
//                     const arma::mat& ddd,
//                     double tol = 1e-6,
//                     int Maxit = 1000,
//                     int grid_points = 50) {
//     double c = a;
//     double fa = beta_func(a, llambdax, y12, ddd);
//     double fb = beta_func(b, llambdax, y12, ddd);
//     double fc = fa;
//     int maxit = Maxit + 1;

//     // First test if we have found a root at an endpoint
//     if (fa == 0.0) {
//         Rcpp::Rcout << "[r_zeroin2] Endpoint solution: beta = " << a << "\n";
//         return {a, 0, 0};
//     }
//     if (fb == 0.0) {
//         Rcpp::Rcout << "[r_zeroin2] Endpoint solution: beta = " << b << "\n";
//         return {b, 0, 0};
//     }

//     // Fallback to grid+Newton if f(a) and f(b) have same sign
//     if (fa * fb > 0.0) {
//         Rcpp::Rcout << "[r_zeroin2] No sign change, falling back to grid+Newton...\n";

//         auto g = [&](double beta) {
//             return beta_func(beta, llambdax, y12, ddd);
//         };
    
//         auto gprime = [&](double beta) {
//             arma::mat log_d = log(ddd);
//             arma::mat term = pow(ddd, beta) % llambdax % log_d;
//             return arma::accu(term % log_d);
//         };

//         // Step 1: coarse grid search
//         double step = (b - a) / (grid_points - 1);
//         double beta_min = a, g_min = std::abs(g(a));
//         for (int i = 1; i < grid_points; i++) {
//             double beta = a + i * step;
//             double val = std::abs(g(beta));
//             if (val < g_min) {
//                 g_min = val;
//                 beta_min = beta;
//             }
//         }

//         // Step 2: Newton refinement
//         double beta = beta_min;
//         for (int it = 0; it < Maxit; it++) {
//             double f = g(beta);
//             double fp = gprime(beta);
//             double update = f / fp;
//             beta -= update;
    
//             if (beta <= a) beta = a + 1e-12;
//             if (beta >= b) beta = b - 1e-12;
    
//             if (std::abs(update) < tol) {
//                 Rcpp::Rcout << "[r_zeroin2] Grid+Newton solution: beta = "
//                             << beta << " (iterations: " << it+1 << ")\n";
//                 arma::vec out(3);
//                 out(0) = beta;
//                 out(1) = 0.0;
//                 out(2) = it+1;
//                 return out;
//             }
//         }

//         Rcpp::warning("Grid+Newton did not converge, returning last iterate");
//         arma::vec out(3);
//         out(0) = beta;
//         out(1) = 0.0;
//         out(2) = -1;
//         return out;
//     }

//     // Otherwise proceed with the Brent loop
//     while (maxit--) { // Main iteration loop
//         double prev_step = b - a;
//         if (std::fabs(fc) < std::fabs(fb)) {
//             a = b; b = c; c = a;
//             fa = fb; fb = fc; fc = fa;
//         }
//         double tol_act = 2 * EPSILON * std::fabs(b) + tol / 2;
//         double new_step = (c - b) / 2;

//         if (std::fabs(new_step) <= tol_act || fb == 0.0) {
//             Rcpp::Rcout << "[r_zeroin2] Brent solution: beta = " << b
//                         << " (iterations: " << Maxit - maxit << ")\n";
//             return {b, std::fabs(c - b), Maxit - maxit};
//         }

//         // Interpolation attempt
//         if (std::fabs(prev_step) >= tol_act && std::fabs(fa) > std::fabs(fb)) {
//             double t1, t2, p, q;
//             double cb = c - b;
//             if (a == c) {
//                 t1 = fb / fa;
//                 p = cb * t1;
//                 q = 1.0 - t1;
//             } else {
//                 q = fa / fc;
//                 t1 = fb / fc;
//                 t2 = fb / fa;
//                 p = t2 * (cb * q * (q - t1) - (b - a) * (t1 - 1.0));
//                 q = (q - 1.0) * (t1 - 1.0) * (t2 - 1.0);
//             }
//             if (p > 0.0) q = -q; else p = -p;
//             if (p < (0.75 * cb * q - std::fabs(tol_act * q) / 2) &&
//                 p < std::fabs(prev_step * q / 2))
//                 new_step = p / q;
//         }

//         if (std::fabs(new_step) < tol_act)
//             new_step = (new_step > 0.0 ? tol_act : -tol_act);

//         a = b; fa = fb;
//         b += new_step; fb = beta_func(b, llambdax, y12, ddd);
//         if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
//             c = a; fc = fa;
//         }
//     }

//     // failed
//     Rcpp::warning("Brent did not converge");
//     return {b, std::fabs(c - b), -1};
// }



// // this is the new code to try out to see if we can pin-point the issue 08/28/2025:
// // [[Rcpp::export]]
// arma::vec r_zeroin2(double a, double b,
//                     const arma::mat& llambdax,
//                     const arma::mat& y12,
//                     const arma::mat& ddd,
//                     double tol=1e-6, int Maxit=1000)
// {
//   // fixed internal params (no API change)
//   const int    max_expand     = 20;
//   const double expand_factor  = 0.05;
//   const double tiny           = 1e-14;

//   auto FAIL = [&](int code, const char* why,
//                   double A, double FA, double B, double FB,
//                   double L, double FL, double R, double FR,
//                   int iters_info)
//   {
//     // One compact, grep-able line to stderr
//     std::fprintf(stderr,
//       "[r_zeroin2 FAIL code=%d] %s | a=%.16g fa=%.16g | b=%.16g fb=%.16g"
//       " | L=%.16g fL=%.16g R=%.16g fR=%.16g | iters=%d\n",
//       code, why, A, FA, B, FB, L, FL, R, FR, iters_info);
//     std::fflush(stderr);
//   };

//   if(!std::isfinite(a) || !std::isfinite(b) || !(a<b)) {
//     FAIL(-4, "bad inputs", a, NA_REAL, b, NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, 0);
//     return arma::vec({NA_REAL, NA_REAL, -4});
//   }

//   double fa = beta_func(a, llambdax, y12, ddd);
//   double fb = beta_func(b, llambdax, y12, ddd);

//   if (bad(fa) || bad(fb)) {
//     FAIL(-3, "NaN/Inf at endpoints", a, fa, b, fb, NA_REAL, NA_REAL, NA_REAL, NA_REAL, 0);
//     return arma::vec({NA_REAL, NA_REAL, -3});
//   }
//   if (near_zero(fa, tiny)) return arma::vec({a, 0.0, 0.0});
//   if (near_zero(fb, tiny)) return arma::vec({b, 0.0, 0.0});

//   // Try to auto-bracket if no sign change; odn't want this because only want to expand to the left to try and find a sign change
//   if (fa*fb > 0.0) {
//     double L = a, R = b, fL = fa, fR = fb;
//     for (int k=0; k<max_expand && fL*fR>0.0; ++k) {
//       double width = (R - L);
//       double newL = std::max(L - (expand_factor-1.0)*width, L - expand_factor*width);
//       double newR = R + expand_factor*width;
//       double fNewL = beta_func(newL, llambdax, y12, ddd);
//       double fNewR = beta_func(newR, llambdax, y12, ddd);
//       if (bad(fNewL) || bad(fNewR)) break;
//       L = newL; R = newR; fL = fNewL; fR = fNewR;
//     }
    
//     std::fprintf(stderr,
//       "[r_zeroin2 TRACE post] a=%.16g fa=%.16g | b=%.16g fb=%.16g | L=%.16g fL=%.16g R=%.16g fR=%.16g\n",
//       a, fa, b, fb, L, fL, R, fR);
//     std::fflush(stderr);
    
//     if (fL*fR > 0.0) {
//       FAIL(-2, "no bracket (after expand)", a, fa, b, fb, L, fL, R, fR, 0);
//       return arma::vec({NA_REAL, NA_REAL, -2});
//     }
//     // adopt expanded bracket if helpful
//     if (fL*fa <= 0.0 && fL*fR < 0.0) { a = L; fa = fL; }
//     if (fR*fb <= 0.0 && fL*fR < 0.0) { b = R; fb = fR; }
//   }

//     // // Try to auto-bracket if no sign change
//     // if (fa * fb > 0.0) {
//     //   // Left-only expansion: keep the right endpoint fixed at its initial value (R=b)
//     //   double L = a, R = b;
//     //   double fL = fa, fR = fb;
    
//     //   for (int k = 0; k < max_expand && fL * fR > 0.0; ++k) {
//     //     double width = (R - L);
//     //     if (width <= 0.0) width = 1.0;                    // safety
//     //     double grow  = (expand_factor - 1.0) * width;     // geometric growth
//     //     double newL  = L - grow;                          // expand only to the left
//     //     double fNewL = beta_func(newL, llambdax, y12, ddd);
//     //     if (bad(fNewL)) break;                            // don't chase NaNs
//     //     L = newL; fL = fNewL;
//     //   }
    
//     //   if (fL * fR > 0.0) {
//     //     // Optional: coarse interior sweep to catch a sign change between L and R
//     //     const int Nprobe = 16;
//     //     double prevx = L, prevf = fL;
//     //     bool found = false;
//     //     for (int i = 1; i <= Nprobe; ++i) {
//     //       double x  = L + (R - L) * (double(i) / double(Nprobe));
//     //       double fx = beta_func(x, llambdax, y12, ddd);
//     //       if (!bad(fx) && prevf * fx <= 0.0) {
//     //         a = prevx; fa = prevf; b = x; fb = fx;        // bracket found inside
//     //         found = true; break;
//     //       }
//     //       prevx = x; prevf = fx;
//     //     }
//     //     if (!found) {
//     //       FAIL(-2, "no bracket (left-only expand)", a, fa, b, fb, L, fL, R, fR, 0);
//     //       return arma::vec({NA_REAL, NA_REAL, -2});
//     //     }
//     //   } else {
//     //     // Success: sign change achieved by left expansion
//     //     a = L; fa = fL; b = R; fb = fR;
//     //   }
//     // }

//   // Brent/zeroin core
//   double c  = a, fc = fa;
//   int maxit = Maxit+1;

//   while (maxit--) {
//     if (std::abs(fc) < std::abs(fb)) {
//       std::swap(a,b); std::swap(b,c);
//       std::swap(fa,fb); std::swap(fb,fc);
//     }

//     const double tol_act = 2.0*macheps()*std::abs(b) + 0.5*tol;
//     double new_step = 0.5*(c - b);

//     if (std::abs(new_step) <= tol_act || fb == 0.0) {
//       return arma::vec({b, std::abs(c-b), double(Maxit - maxit)});
//     }

//     double prev_step = b - a;

//     if (std::abs(prev_step) >= tol_act && std::abs(fa) > std::abs(fb)) {
//       double p, q;
//       if (a == c) {
//         double denom = (fa - fb);
//         if (denom == 0.0) denom = (denom >= 0 ? tiny : -tiny);
//         p = (b - a) * fb / denom; q = 1.0;
//       } else {
//         double q1 = fa/fc, t1 = fb/fc, t2 = fb/fa;
//         p = t2 * ((c-b)*q1*(q1 - t1) - (b - a)*(t1 - 1.0));
//         q = (q1 - 1.0)*(t1 - 1.0)*(t2 - 1.0);
//       }
//       if (p > 0) q = -q; else p = -p;
//       if (p < (0.75*(c-b)*q - 0.5*std::abs(tol_act*q)) &&
//           p < 0.5*std::abs(prev_step*q)) {
//         new_step = p/q;
//       }
//     }

//     if (std::abs(new_step) < tol_act) new_step = (new_step > 0 ? tol_act : -tol_act);

//     a = b;  fa = fb;
//     b += new_step;
//     fb = beta_func(b, llambdax, y12, ddd);
//     if (bad(fb)) {
//       FAIL(-3, "NaN/Inf during iter", a, fa, b, fb, NA_REAL, NA_REAL, NA_REAL, NA_REAL,
//           Maxit - maxit);
//       return arma::vec({NA_REAL, NA_REAL, -3});
//     }

//     if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) { c = a; fc = fa; }
//   }

//   FAIL(-1, "max iterations", a, fa, b, fb, NA_REAL, NA_REAL, NA_REAL, NA_REAL, Maxit);
//   return arma::vec({b, std::abs(c-b), -1});
// }



///////////////////////////////////
///////Project onto the cube
///////////////////////////////////
// [[Rcpp::export]]
arma::vec project_func(const arma::vec& a,const arma::mat& bounds) {
    arma::vec a_copy(a.n_elem);
    for (uint i=0; i<bounds.n_rows; i++) {
        if ((a(i) > bounds(i,0)) && (a(i)<bounds(i,1))) {
            a_copy(i) = a(i);
        } else {
            if (fabs(a(i) - bounds(i, 0)) < fabs(a(i) - bounds(i, 1)) ) {
                a_copy(i) = bounds(i,0);
            } else {
                a_copy(i) = bounds(i,1);
            }

        }
    }
    return a_copy;
}

//return_grad(y_k, S1, S2, llambdax, y12)
// res = grad(nll, x0, u, v, llambdax, y12, 0.001, true, npt);


///////////////////////////////////
///////objective function: negative Log Likelihood Function
///////////////////////////////////
// [[Rcpp::export]]
double nll(const arma::vec& phi,const  arma::mat& u,const  arma::mat& v,const  arma::mat& llambdax,const  arma::mat& y12) {
    arma::mat ddd = dist_trans(phi, u, v);
    arma::vec beta_sol = r_zeroin2(-2, 0, llambdax, y12, ddd);
    double beta = beta_sol(0);
    double like = arma::accu(pow(ddd,beta) % llambdax - beta * y12 % log(ddd));
    return like;
}


// gradient, analytical
// [[Rcpp::export]]
arma::vec grad(const arma::vec& phi, const arma::mat& u, const arma::mat& v,
                const arma::mat& llambdax, const arma::mat& y12){
    arma::mat ddd = dist_trans(phi, u, v);
    arma::vec beta_sol = r_zeroin2(-2, 0, llambdax, y12, ddd);
    double beta = beta_sol(0);
    arma::mat part1 = pow(ddd,beta) % llambdax;
    arma::mat p1(3,3), p2(3,3), p3(3,3);
    p1(1,2) = -1;
    p1(2,1) = 1;
    p2(0,2) = 1;
    p2(2,0) = -1;
    // p3(0,1) = -1;
    // p3(1,0) = 1;
    p3(0,1) = 1;
    p3(1,0) = -1;
    arma::vec res(3);
    // arma::mat result = -2 * x * Rz(phi(2)) * Ry(phi(1)) * Rx(phi(0)) * trans(y);
    arma::mat temp = -2 * u * Rz(phi(2)) * Ry(phi(1)) * Rx(phi(0)) * p1 * arma::trans(v);
    temp = temp / ddd;
    res(0) = beta *accu(temp% part1 - temp%y12);
    temp = -2 * u * Rz(phi(2)) * Ry(phi(1)) * p2 * Rx(phi(0))* arma::trans(v);
    temp = temp / ddd;
    res(1) = beta *accu(temp% part1 - temp%y12);
    temp = -2 * u * Rz(phi(2)) * p3 *Ry(phi(1)) *  Rx(phi(0))* arma::trans(v);
    temp = temp / ddd;
    res(2) = beta *accu(temp% part1 - temp%y12);

    return res;
}

// // [[Rcpp::export]]
// double beta_func(const arma::mat& d1, const arma::mat& d2, const arma::mat& log_d)) {
//     double y = arma::accu(log_d % d1 - d2);
//     return y;
// }

// // pow(ddd, beta) % llambdax
// // [[Rcpp::export]]
// inline arma::mat d_beta_lambda(const double& beta, const arma::mat& llambdax, const arma::mat& ddd) {
//     arma::mat res = pow(ddd,beta) %llambdax;
//     return res;
// }

// // y12 % log_d
// // [[Rcpp::export]]
// inline arma::mat y_d(const arma::mat& y12, const arma::mat& log_d) {
//     arma::mat res=y12 % log_d;
//     return res;
// }




// double nll1(const arma::vec& phi,const  arma::mat& u,const  arma::mat& v,const  arma::mat& llambdax,const  arma::mat& y12) {
//     arma::mat ddd = dist_trans(phi, u, v);
//     arma::vec beta_sol = r_zeroin2(-2, 0, llambdax, y12, ddd);
//     double beta = beta_sol(0);
//     double like = arma::accu(pow(ddd,beta) % llambdax - beta * y12 % log(ddd));
//     return like;
// }

// 1 iteration
// [[Rcpp::export]]
arma::vec APG(const arma::mat& y12,
              const arma::mat& S1,
              const arma::mat& S2,
              const arma::mat& llambdax,
              int T=500,
              double tol=0.0000001,
              double rel_tol=0.01)
{
    double n = y12.n_rows * y12.n_cols;
    arma::vec a1 = rcpp_gen_angle();
    // arma::vec a1 = {0.3343907, 4.3702400, 2.3869913};
    arma::mat bounds = {{0, 2*M_PI},
                        {0, 2*M_PI},
                        {0, M_PI} };
    double diff = 10, x_diff = 10;
    arma::vec t(T), q(T), c_value(T), likelihood(T);
    t(0) = 0;
    t(1) = 1;
    q(0) = 1;
    q(1) = 1;
    double a_x = 0.1/n;
    double a_y = 0.1/n;
    arma::mat x(3,T);
    // vec to row vec
    x.col(0) = a1;
    x.col(1) = a1;
    int s=1;
    double theta = 0.01, ita = 0.1;
    arma::vec z_cur(3), z_next(3), v_next(3), scaled_y_diff(3), temp(3), y_k(3);
    z_cur = a1;
    double obj_z_next, diff1;
    c_value(0) = nll(a1, S1, S2, llambdax, y12);
    c_value(1) = nll(a1, S1, S2, llambdax, y12);

    while ((s<T-1) & ((diff > tol) | (x_diff > rel_tol))) {
        // everything is column based
        y_k = x.col(s) - t(s-1) / t(s) * (z_cur - x.col(s)) + (t(s-1) - 1) / t(s) * (x.col(s) - x.col(s-1));
        scaled_y_diff = y_k - a_y * grad(y_k, S1, S2, llambdax, y12);
        // z_next <- project_function(temp, bounds)
        z_next = project_func(scaled_y_diff, bounds);
        // obj_func(a = z_next, u = S1, v = S2, llambdax = llambdax, y12 = y12)
        obj_z_next = nll(z_next, S1, S2, llambdax, y12);
        diff1 = obj_z_next - c_value(s) + theta * arma::sum(arma::square(z_next-y_k));
        if (diff1 < 0) {
            x.col(s+1) = z_next;
        } else {
            temp = x.col(s) - a_x * grad(x.col(s), S1, S2, llambdax, y12);
            v_next = project_func(temp, bounds);
            if (obj_z_next < nll(v_next, S1, S2, llambdax, y12)) {
                x.col(s+1) = z_next;
            } else {
                x.col(s+1) = v_next;
            }
        }
        t(s+1) = (sqrt(arma::as_scalar(4*t(s)*t(s) + 1)) + 1)/2;
        q(s+1) = ita * q(s) + 1;
        // nll for x(s+1)
        likelihood(s+1) = nll(x.col(s+1), S1, S2, llambdax, y12);
        c_value(s+1) = (ita * q(s) * c_value(s) + likelihood(s+1))/q(s+1);
        z_cur = z_next;
        diff = fabs(likelihood(s+1) - likelihood(s));
        temp = arma::abs((x.col(s+1) - x.col(s))/ x.col(s));
        // fill NA with 0
        temp.transform( [](double val) { return (std::isnan(val) ? double(0) : val); } );
        x_diff = temp.max();
        s++;
    }

    arma::vec res(6);
    arma::mat ddd = dist_trans(x.col(s), S1, S2);
    arma::vec beta_sol = r_zeroin2(-2, 0, llambdax, y12, ddd);

    for (int i=0; i<3; i++) {
        res(i) = x(i,s);
    }
    res(3) = likelihood(s);
    res(4) = beta_sol(0);
    res(5) = s;
    return res;
}



// [[Rcpp::export]]
arma::vec minimizer(const uint Nrep,
                    const arma::mat& y12,
                    const arma::mat& S1,
                    const arma::mat& S2,
                    const arma::mat& llambdax,
                    int T=500,
                    double tol=0.0000001,
                    double rel_tol=0.01,
                    int threads=1)
{
    #ifdef _OPENMP
    if ( threads > 0 )
        omp_set_num_threads( threads );
    // REprintf("Number of threads=%i\n", omp_get_max_threads());
    #endif

    arma::vec a1={0,0,0};
    // both arma and NumericVector are column major
    arma::mat res(6,Nrep);
    #pragma omp parallel for schedule(dynamic)
    for (uint i=0; i<Nrep; i++) {
        res.col(i) = APG(y12, S1, S2, llambdax, T, tol, rel_tol);
    //    printf("Thread %d is running number %d\n", omp_get_thread_num(), i);
    }

    arma::vec result = res.col(res.row(3).index_min());
    // check how many times (out of Nrep) the minimum have been reached
    arma::uvec min_reached = arma::find(arma::abs(res.row(3) - res.row(3).min()) < 0.01);
    // grow by 1 to record length(min_reached)
    result.resize(result.n_elem+1);
    result.back() = min_reached.n_elem;
    return result;
}
