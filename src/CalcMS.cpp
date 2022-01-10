#include<Rcpp.h> //On va utiliser les fonctions et classes de l'en-tete du fichier Rcpp.h
#include<omp.h> //Pour paralléliser facilement du code C++
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
double moy(NumericVector x) {
 double total = 0;
 int n = x.size();
 for (int i=0; i<n; i++) {
   total += x[i];
 }
 return total/n;
}


// [[Rcpp::export]]
double ect(NumericVector x) {
 double total = 0;
 double m = 0;
 int n = x.size();
 for (int i=0; i<n; i++) {
   total += pow(x[i],2);
   m += x[i];
 }
 return sqrt(total/n - pow(m/n,2));
}


// [[Rcpp::export]]
double moy2(NumericVector x) {
 double total = 0;
 int n = x.size();
 #pragma omp parallel for
 for (int i=0; i<n; i++) {
   total += x[i];
 }
 return total/n;
}


// [[Rcpp::export]]
double ect2(NumericVector x) {
 double total = 0;
 double m = 0;
 int n = x.size();
 #pragma omp parallel for
 for (int i=0; i<n; i++) {
   total += pow(x[i],2);
   m += x[i];
 }
 return sqrt(total/n - pow(m/n,2));
}
