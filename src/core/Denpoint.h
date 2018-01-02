#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class Denpoint {
public:
	vector<double> point;
	bool covered;
	int timestamp;

	Denpoint (vector<double> point,int timestamp) {
		this->point = point;
		this->covered = false;
		this->timestamp = timestamp;
	}

	public:void set_covered () {
		this->covered = !this->covered;
	}
};