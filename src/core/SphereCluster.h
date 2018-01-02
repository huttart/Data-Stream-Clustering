#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class Sphere {
public:
	int N;
	double r;
	vector<double> center;
	vector<double> LS;

	Sphere (double r,vector<double> center) {
		this->r = r;
		this->center = center;
		this->LS = center;
		N = 1;
	}
	
	public:void insert (vector<double> point) {
		N++;
		for (int i = 0 ; i < point.size() ; i++) {
			LS[i] += point[i];
			center[i] = LS[i]/N;
		}
	}
};