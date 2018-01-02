#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

double dist(NumericVector point1,NumericVector point2) {
	int size = point1.size();
	double distance = 0;
	double tmp = 0;
	for (int i = 0 ; i < size ; i++) {
		tmp = point1[i] - point2[i];
		distance += tmp*tmp;
	}

	return sqrt(distance);
}

// [[Rcpp::export]]
int reassign (NumericVector x,List center) {
	int size = center.size();
	double mindist = 1000000;
	
	int index = 0;
	for (int i = 0 ; i < size ; i++) {
		NumericVector tmp = center[i];
		double distance = dist(x,tmp);
		// cout << distance << endl;
		if (distance < mindist) {
			mindist = distance;
			index = tmp[x.size()];
		}
	}
	return index;
}