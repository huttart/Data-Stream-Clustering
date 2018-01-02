#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class MicroCluster {
private:
	int	N;
	int lastEdit;		
    int	creationTime;					
    double	lambda;
    vector<double>	LS; 				
    vector<double>	SS; 
public:
	bool covered; // this variable is used in finalCluster function
    vector<vector<double> > datapoints;
	MicroCluster (vector<double> point,double lambda,int timastamp)
	{
		N = 1;
		lastEdit = timastamp;
		creationTime = timastamp;
		this->lambda = lambda;
		int size = point.size();

		for (int i = 0 ; i < size ; i++) {
			LS.push_back(point[i]);
			SS.push_back(point[i]*point[i]);
		}
		this->covered = false;
		point.push_back(timastamp);
		datapoints.push_back(point);
	}	

	public:void insert_point(vector<double> point,int timastamp)
	{
		int size = point.size();
		// cout << "size = " << size << endl;
		for(int i = 0 ; i < size ; i++)
		{
			// cout << point[i] << " ";
			LS[i] += point[i];
			SS[i] += point[i]*point[i];
		}

		N += 1;
		lastEdit = timastamp;
		point.push_back(timastamp);
		datapoints.push_back(point);
		// cout << endl;
	}	

	public:double getweight(int currentTimestamp)
	{
		int dt = currentTimestamp - lastEdit;
		return N*(pow(2,(-1)*lambda*dt));
	}

	public:vector<double> getcenter(int currentTimestamp)
	{
		double weight;
		int length = LS.size();
		int dt = currentTimestamp - lastEdit;
		vector<double> center;

		weight = this->getweight(currentTimestamp);


		for(int i = 0 ; i < length ; i++)
		{
			center.push_back(LS[i]*pow(2,(-1)*lambda*dt));
			center[i] = center[i]/weight;
		}

		return center;
	}

	public:vector<double> getCF1(int currentTimestamp)
	{
		int dt = currentTimestamp - lastEdit;
		int length = LS.size();
		vector<double> CF1;

		for(int i = 0 ; i < length ; i++)
		{
			CF1.push_back(LS[i]*pow(2,(-1)*lambda*dt));
		}

		return CF1;

	}

	public:vector<double> getCF2(int currentTimestamp)
	{
		int dt = currentTimestamp - lastEdit;
		int length = SS.size();
		vector<double> CF2;

		for(int i = 0 ; i < length ; i++)
		{
			CF2.push_back(SS[i]*pow(2,(-1)*lambda*dt));
		}

		return CF2;

	}

	public:double getradius(int currentTimestamp) 
	{
		double weight;

		int length = SS.size();
		int dt = currentTimestamp - lastEdit;

		vector<double> cf2;
		vector<double> cf1;

		weight = this->getweight(currentTimestamp);
		cf1 = this->getCF1(currentTimestamp);
		cf2 = this->getCF2(currentTimestamp);

		double max = 0;

	    for (int i = 0; i < length; i++) {
	        double x1 = cf2[i] / weight;
	        double x2 = pow(cf1[i] / weight, 2);
	        double tmp = x1 - x2;

	        if(tmp < 0){
	        	continue;
	        }

	        if (sqrt(tmp) > max) {
	            max = sqrt(tmp);
	        }
	    }
	    // cout << "radius = " << max << endl;
	    return max;

	}

	public:int getCreationTime()
	{
		return this->creationTime;
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


