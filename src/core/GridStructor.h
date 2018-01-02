#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class GS {
public:
	int n; // number of points in grid
	int t; // last update time
	double w; // grid weight
	double lambda;
	double max_distance;
	vector<double> LS; // sum of data points
	vector<vector<double> > data_points;
	GS () {}
	GS (int timestamp,double lambda,vector<double> input) {
		this->n = 1;
		this->t = timestamp;
		this->w = 1;
		this->lambda = lambda;
		this->LS = input;
		this->data_points.push_back(input);
	}

	public:void update_vector (int timestamp,vector<double> input,double lambda) {
		
		this->n += 1;
		this->w *= pow(2,-this->lambda*(timestamp - this->t));
		this->w++;
		this->t = timestamp;
		vector<double> temp = input;
		temp.push_back(timestamp);
		this->data_points.push_back(temp);

		vector<double> center;
		for(int i = 0 ; i < input.size() ; i++) {
			LS[i] += input[i];
			center.push_back(LS[i]/n);
		}
		
		double max = 1000000;

		for (int j = 0 ; j < data_points.size() ; j++) {
			double distance = 0;
			for (int i = 0 ; i < center.size() ; i++) {
				double tmp = data_points[j][i] - center[i];
				distance += tmp*tmp;
			}
			if (sqrt(distance) > max_distance) {
				max_distance = sqrt(distance);
			}
		}


	}
};