#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class MiniCluster {
private:
	int	N;
	int lastEdit;		
    double	lambda;
    double w;
    double r;
    double max_distance;
    vector<double> LS; 				
public:
	bool covered; // this variable is used in finalCluster function
	MiniCluster (vector<double> LS,double lambda,int timastamp,int n,double max_distance,double w,vector<vector<double> > data_points)
	{
		this->N = n;
		this->lastEdit = timastamp;
		this->lambda = lambda;
		this->LS = LS;
		this->max_distance = max_distance;
		this->w = w;

		vector<double> center = getcenter();
		double sum = 0 ;
		for (int i = 0 ; i < data_points.size() ; i++) {
			double dis = 0;
			for (int j = 0 ; j < center.size() ; j++) {
				double tmp = center[j] - data_points[i][j];
				dis += tmp*tmp;
			}
			sum += sqrt(dis);
		}

		this->r = sum/data_points.size() ;
	}	

	public:void insert_point(vector<double> point,int timastamp)
	{
		int size = point.size();
		for(int i = 0 ; i < size ; i++) {
			LS[i] += point[i];
		}

		N += 1;
		lastEdit = timastamp;

		w = w*pow(2,-lambda*(timastamp - lastEdit)) + 1;
		vector<double> center = getcenter();
		double dis = 0;
		for(int i = 0 ; i < size ; i++) {
			double tmp = center[i] - point[i];
			dis += tmp*tmp;
		}

		this->r += sqrt(dis)/N ;

	}	

	public:double getweight(int currentTimestamp)
	{
		int dt = currentTimestamp - lastEdit;
		// cout << "dt = " << lastEdit << endl;
		// cout << "real_w = " << w << endl;
		// cout << "w = " << w*(pow(2,(-1)*lambda*dt)) << endl;

		// return w*(pow(2,-lambda*dt))+1;
		return w;
	}

	public:vector<double> getcenter()
	{
		int length = LS.size();
		vector<double> center;

		for(int i = 0 ; i < length ; i++) {
			center.push_back(LS[i]/N);
		}

		return center;
	}
	
	public:double getradius() 
	{
		return this->r ;
	}

	public:double get_maxdis() {
		return this->max_distance;
	}

	public:int get_num() {
		return this->N ;
	}


	public:void print()
	{
		cout << "LS: { ";
		for(int i = 0 ; i < LS.size() ; i++)
		{
			cout << LS[i] << " ";
		}

		cout << "}\n";
	}
   
};


