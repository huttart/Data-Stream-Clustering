#include"core/Clustering.h"
#include"core/MicroCluster.h"
#include"core/Denpoint.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


class Denstream {
private:
	int minpts;
	int timestamp;
	int init_points;
	int tp;
	int numberOfCluster;
	int *cl;

	double lambda;
	double mu;
	double epsilon;
	double beta;
	double offline;

	bool initialized;

	Clustering<MicroCluster> p_micro_cluster;
	Clustering<MicroCluster> o_micro_cluster;

	vector<Denpoint> initBuffer;
	vector<vector<double> > missing_points;

public:

	Denstream(double epsilon,double mu,int init_points,double lambda,double beta,int offline) {
		this->epsilon 		= 	epsilon;
		this->mu 			= 	mu;
		this->init_points	=	init_points;
		this->minpts		=	mu;
		this->lambda		=	lambda;
		this->beta 			=	beta;
		this->offline		=	offline;

		if (beta * mu > 1) {
			tp = round(1 / lambda * log(((beta * mu)/(beta * mu - 1)))) + 1;
		} else {
			tp = 1;
		}

		cout << "tp = " << tp << endl;

		timestamp = 0;
		initialized = false;
		numberOfCluster = 1;
	}

	private:double dist(vector<double> point1,vector<double> point2) {
		int size = point1.size();
		double distance = 0;
		double tmp = 0;
		for (int i = 0 ; i < size ; i++) {
			tmp = point1[i] - point2[i];
			distance += tmp*tmp;
		}

		return sqrt(distance);
	}

	private:int nearestMC(vector<double> input,Clustering<MicroCluster> x) {
		int min = 0;
		double distance = 1000000;
		int size = x.sizes();

		for (int i = 0 ; i < size ; i++) {
			MicroCluster mc = x.get(i);
			double dis = dist(mc.getcenter(timestamp),input) - mc.getradius(timestamp);

			if (dis < distance) {
				distance = dis;
				min = i;
			}
		}

		return min;	

	}

	private:bool try_merge(vector<double> input,Clustering<MicroCluster> *x,int index) {
		bool merge = false;
		MicroCluster mc = x->get(index);
		
		mc.insert_point(input,timestamp);

		if (mc.getradius(timestamp) <= epsilon) {
			x->getreal(index)->insert_point(input,timestamp);
			merge = true;
		}

		return merge;
	}

	private:void remove() {
		int size = p_micro_cluster.sizes();
		vector<int> remove_list;

		for (int i = 0 ; i < size ; i++) {
			MicroCluster mc = p_micro_cluster.get(i);
			if (mc.getweight(timestamp) < (beta * mu)) {
				remove_list.push_back(i);
			}
		}

		for (int i = 0 ; i < remove_list.size() ; i++) {
			p_micro_cluster.remove(remove_list[i]);
		}

		size = o_micro_cluster.sizes();
		vector<int> remove_list2;

		for (int i = 0 ; i < size ; i++) {
			MicroCluster mc = o_micro_cluster.get(i);
			int t0 = mc.getCreationTime();
            double xsi1 = pow(2, (-lambda * (timestamp - t0 + tp))) - 1;
            double xsi2 = pow(2, -lambda * tp) - 1;
            double xsi = xsi1 / xsi2;

			if (mc.getweight(timestamp) < xsi) {
				remove_list2.push_back(i); 
			}
		}

		for (int i = 0 ; i < remove_list2.size() ; i++) {
			o_micro_cluster.remove(remove_list2[i]);
		}

	}

	private:void initialDBScan() {
		for (int i = 0 ; i < initBuffer.size() ; i++) {
			Denpoint *point = &initBuffer[i];
			if (!point->covered) {
				point->covered = true;
				vector<int> neighbourhood = getNeighbourhoodIDs(point, &initBuffer, epsilon);
				if (neighbourhood.size() > minpts) {
					MicroCluster mc(point->point, lambda,  point->timestamp);
					expandCluster(&mc, &initBuffer, neighbourhood);
					p_micro_cluster.add(mc);
				} else {
					point->covered = false;
				}
			}
		}
	}

	private:void  expandCluster(MicroCluster *mc, vector<Denpoint> *points, vector<int> neighbourhood) {
		for (int i = 0 ; i < neighbourhood.size() ; i++) {
			Denpoint *point = &(*points)[neighbourhood[i]];
			if (!point->covered) {
				point->covered = true;
				mc->insert_point(point->point, point->timestamp);
				vector<int> neighbourhood2 = getNeighbourhoodIDs(point, points, epsilon);
				if (neighbourhood2.size() > minpts) {
					expandCluster(mc, points, neighbourhood2);
				}
			}
		}
	}

	private:vector<int> getNeighbourhoodIDs (Denpoint *point, vector<Denpoint> *points, double eps) {
		vector<int> neighbourhood;
		for (int i = 0 ; i < points->size() ; i++) {
			Denpoint point2 = (*points)[i];
			if (!point2.covered) {
				double distance = dist(point->point, point2.point);
				if (distance < eps) {
					neighbourhood.push_back(i);
				}
			}
		}

		return neighbourhood;
	}

	private:vector<double> convert_input(NumericVector x) {
		int size = x.size();
		vector<double> tmp;
		for (int i = 0 ; i < size ; i++) {
			tmp.push_back(x[i]);
		}
		return tmp;
	}

	private:void missingPoints () {
		for (int i = 0 ; i < initBuffer.size() ; i++) {
			Denpoint point = initBuffer[i];
			if(!point.covered){
				vector<double> tmp = point.point;
				tmp.push_back(point.timestamp);
				missing_points.push_back(tmp);
			}
		}
	}

	public:void trainOnInstance(NumericVector point) {
		timestamp++;
		// cout << "timestamp = " << timestamp << endl;
		vector<double> input = convert_input(point);
		bool merge = false;
		int nearest = 0;

		if (!initialized) {
			initBuffer.push_back(Denpoint(input,timestamp));
            if (initBuffer.size() >= init_points) {
                initialDBScan();
                initialized = true;
                missingPoints();
            }
		} else {
			if(p_micro_cluster.sizes() > 0) {
				nearest = nearestMC(input,p_micro_cluster);
				merge = try_merge(input,&p_micro_cluster,nearest);

			}
			if (!merge && o_micro_cluster.sizes() > 0) {
				nearest = nearestMC(input,o_micro_cluster);
				merge = try_merge(input,&o_micro_cluster,nearest);



				if (merge) {
					MicroCluster mc = o_micro_cluster.get(nearest);
					if (mc.getweight(timestamp) > (beta * mu)) {
						// cout << "weight: " << mc.getweight(timestamp) << endl;
						p_micro_cluster.add(mc);
						o_micro_cluster.remove(nearest);
					}
				}
			}
			if (!merge) {	
				MicroCluster mc(input, lambda, timestamp);
				o_micro_cluster.add(mc);
			}
			if (timestamp%tp == 0) {
				remove();
			}
		}

	}

	public:List final()
	{
		return finalCluster(vector<double>(),epsilon*offline,true); 
	}

	public:List finalCluster(const vector<double> center = vector<double>(),double eps = 0,bool starting = false) {
		// cout << "epsilon = " << eps << endl;
		if (starting && p_micro_cluster.sizes() > 0) {
			cl = new int[p_micro_cluster.sizes()];
			for (int i = 0 ; i < p_micro_cluster.sizes() ; i++) {
				cl[i] = 0;
			}
			for (int i = 0 ; i < p_micro_cluster.sizes() ; i++) {
				MicroCluster mc = p_micro_cluster.get(i);
				if (!cl[i]) {
					numberOfCluster += 1;
					cl[i] = numberOfCluster;
					finalCluster(mc.getcenter(timestamp),mc.getradius(timestamp));
				}
			}
		} else {
			for (int i = 0 ; i < p_micro_cluster.sizes() ; i++) {
				MicroCluster mc = p_micro_cluster.get(i);
				if (!cl[i]) {
					vector<double> center2 = mc.getcenter(timestamp);
					double distance = dist(center,center2);
					if (distance <= eps+mc.getradius(timestamp)) {
						cl[i] = numberOfCluster;
						finalCluster(center2,mc.getradius(timestamp));
					}
				}
			}
		}
		List microCluster(p_micro_cluster.sizes());
		int n = 0;
		for (int i = 0 ; i < p_micro_cluster.sizes() ; i++) {
			MicroCluster mc = p_micro_cluster.get(i);
			vector<double> tmp = mc.getcenter(timestamp);
			tmp.push_back(cl[i]);
			tmp.push_back(mc.getradius(timestamp));
			tmp.push_back(0);
			microCluster[n++] = tmp;
		}
	
		return microCluster;
	}

	private:List centroids_of_cluster(int num_of_cluster,int *cluster,Clustering<MicroCluster> pmc) {
		List centroids(num_of_cluster);

		cout << "num = " << num_of_cluster << endl;
		for (int i = 0 ; i < num_of_cluster ; i++) {
			vector<double> tmp;
			int n = 0;
			for (int j = 0 ; j < pmc.sizes() ; j++) {
				if(cluster[j] == i+1) {
					MicroCluster mc = pmc.get(j);
					vector<double> center = mc.getcenter(timestamp);
					if(tmp.size() == 0) {
						for (int k = 0 ; k < center.size() ; k++) {
							tmp.push_back(center[k]);
						}
					}
					for (int k = 0 ; k < center.size() ; k++) {
						tmp[k] += center[k];
					}
					n++;
				}
			}
			for (int k = 0 ; k < tmp.size() ; k++) {
				tmp[k] /= n;
			}
			centroids[i] = tmp;
		}

		return centroids;
	}

	public:List end() {
		List x = final();
		cout << "\nmicro-cluster: " << p_micro_cluster.sizes() << endl;
		cout << "marcro-cluster: " << numberOfCluster-1 << endl;
		// cout << "timestamp: " << timestamp << endl;

		return x;
	}

	public:List reassign () {
		List result(timestamp);
		int n = 0;
		for (int i = 0 ; i < p_micro_cluster.sizes() ; i++) {
			MicroCluster mc = p_micro_cluster.get(i);
			for (int j = 0 ; j < mc.datapoints.size() ; j++) {
				vector<double> tmp;
				tmp.push_back(mc.datapoints[j][mc.datapoints[j].size()-1]);
				tmp.push_back(cl[i]);
				result[n++] = tmp;
			}
			
		}
		for (int i = 0 ; i < o_micro_cluster.sizes() ; i++) {
			MicroCluster mc = o_micro_cluster.get(i);
			for (int j = 0 ; j < mc.datapoints.size() ; j++) {
				vector<double> tmp = mc.datapoints[j];
				tmp.erase(tmp.begin() + tmp.size()-1);
				double mindist = 1000000;
				int cluster = 0;
				for (int k = 0 ; k < p_micro_cluster.sizes() ; k++) {
					double distance = dist(p_micro_cluster.get(k).getcenter(timestamp),tmp);
					if (distance < mindist) {
						mindist = distance;
						cluster = 999999;
					}
				}
				vector<double> tmp2;
				tmp2.push_back(mc.datapoints[j][mc.datapoints[j].size()-1]);
				tmp2.push_back(cluster);
				result[n++] = tmp2;
			}
			
		}

		for (int i = 0 ; i < missing_points.size() ; i++) {
			vector<double> tmp;
			tmp.push_back(missing_points[i][missing_points[i].size() - 1]);
			tmp.push_back(999999);
			result[n++] = tmp;
		}

	
		return result;
	}

	public:List reassign2 () {
		List result(timestamp);
		int n = 0;
		for (int i = 0 ; i < p_micro_cluster.sizes() ; i++) {
			MicroCluster mc = p_micro_cluster.get(i);
			for (int j = 0 ; j < mc.datapoints.size() ; j++) {
				vector<double> tmp;
				tmp.push_back(mc.datapoints[j][mc.datapoints[j].size()-1]);
				tmp.push_back(cl[i]);
				result[n++] = tmp;
			}
			
		}
		for (int i = 0 ; i < o_micro_cluster.sizes() ; i++) {
			MicroCluster mc = o_micro_cluster.get(i);
			for (int j = 0 ; j < mc.datapoints.size() ; j++) {
				vector<double> tmp = mc.datapoints[j];
				tmp.erase(tmp.begin() + tmp.size()-1);
				double mindist = 1000000;
				int cluster = 0;
				for (int k = 0 ; k < p_micro_cluster.sizes() ; k++) {
					double distance = dist(p_micro_cluster.get(k).getcenter(timestamp),tmp);
					if (distance < mindist) {
						mindist = distance;
						cluster = cl[k];
					}
				}
				vector<double> tmp2;
				tmp2.push_back(mc.datapoints[j][mc.datapoints[j].size()-1]);
				tmp2.push_back(cluster);
				result[n++] = tmp2;
			}
			
		}

		for (int i = 0 ; i < missing_points.size() ; i++) {
			vector<double> tmp;
			tmp.push_back(missing_points[i][missing_points[i].size() - 1]);
			tmp.push_back(999999);
			result[n++] = tmp;
		}

	
		return result;
	}





};


RCPP_EXPOSED_CLASS(Denstream);

RCPP_MODULE(yada) {
	using namespace Rcpp;

  	class_<Denstream>("Denstream")
    	.constructor< double, double, int, double, double, int  >()
    	.method("trainOnInstance", &Denstream::trainOnInstance)
    	.method("end", &Denstream::end)
    	.method("final", &Denstream::final)
    	.method("reassign", &Denstream::reassign)
    	.method("reassign2", &Denstream::reassign2)
	;
};




