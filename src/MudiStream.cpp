#include"core/CoreMiniCluster.h"
#include"core/GridStructor.h"
#include"core/QuickSort.h"	
#include"core/Clustering.h"	
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

class Mudi {
private:
	int N;
	int tp;
	int timestamp;
	int minpts;
	double alpha;
	double lambda;
	double gridSize;
	int *cl;
	map<vector<int>, GS> grid_list;
	Clustering<MiniCluster> cmc;
public:
	Mudi (int n,double lambda,double alpha,double gridSize,int minpts) {
		this->N = n;
		this->lambda = lambda;
		this->alpha = alpha;
		this->gridSize = gridSize;
		this->minpts = minpts;	
		if ( (alpha - (n * (1 - pow(2,-lambda) ) ) ) > 0) {
			this->tp = round(1 / lambda * log( alpha/(alpha - (n * (1 - pow(2,-lambda) ) ) ) ) );
		} else {
			this->tp = 1;
		}
		cout << "tp = " << tp << endl;
		cout << "minweight = " << (alpha/(N*(1-pow(2,-lambda)))) << endl;
		this->timestamp = 0;
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

	private:int nearestMC(vector<double> input,Clustering<MiniCluster> x) {
		int min = 0;
		double distance = 1000000;
		int size = x.sizes();
		for (int i = 0 ; i < size ; i++) {
			MiniCluster mc = x.get(i);
			double dis = dist(mc.getcenter(),input) - mc.getradius();

			if (dis < distance) {
				distance = dis;
				min = i;
			}
		}
		return min;	
	}

	public:void trainOnInstance(NumericVector point) {
		vector<double> input = as<vector<double> >(point);
		timestamp++;
		int nearest;
		bool merge = false;
		if (cmc.sizes() > 0) {
			nearest = nearestMC(input,cmc);
			MiniCluster mc = cmc.get(nearest);
			if (dist(input,mc.getcenter()) <= mc.get_maxdis()) {
				cmc.getreal(nearest)->insert_point(input,timestamp);
				merge = true;
			}
		}

		if (!merge) {
			vector<int> grid;
			GS *mc;
			for (int i = 0 ; i < input.size() ; i++) {
				grid.push_back(floor(input[i]/gridSize));
			}
			if (grid_list.insert(make_pair(grid,GS(timestamp,lambda,input))).second == false) {
				mc = &grid_list[grid];
				mc->update_vector(timestamp,input,lambda);
			} else {
				mc = &grid_list[grid];
			}
			if (mc->n >= 1 && mc->w >= (alpha/(N*(1-pow(2,-lambda)))) ) {
				cmc.add(MiniCluster(mc->LS, lambda, timestamp, mc->n, mc->max_distance, mc->w, mc->data_points));
				grid_list.erase(grid);
			}
		}

		if (timestamp%tp == 0) {
			map<vector<int>, GS>::iterator tmp = grid_list.begin();
			while (tmp != grid_list.end()) {
				tmp->second.w *= pow(2,-lambda*(timestamp - tmp->second.t));
				double xs1 = alpha*(1-pow(2,-lambda*(timestamp - tmp->second.t + 1)));
				double xs2 = N*(1-pow(2,-lambda*tmp->second.t));
				double xs = xs1/xs2;
				if (tmp->second.w < xs) {
					grid_list.erase(tmp++);
				} else {
					++tmp;
				}
			}
			vector<int> remove_list;

			for (int i = 0 ; i < cmc.sizes() ; i++) {
				MiniCluster mc = cmc.get(i);
				if (mc.getweight(timestamp) < (alpha/(N*(1-pow(2,-lambda)))) ) {
					remove_list.push_back(i);
				}
			}
			for (int i = 0 ; i < remove_list.size() ; i++) {
				cmc.remove(remove_list[i]);
			}
		}

	}

	private:vector<int> cmcNeighbourhood (Clustering<MiniCluster> x,MiniCluster mc) {
		vector<int> neighbour;
		for (int i = 0 ; i < x.sizes() ; i++) {
			MiniCluster tmp = x.get(i);
			vector<double> cen_tmp = tmp.getcenter();
			vector<double> center = mc.getcenter();
			for (int j = 0 ; j < center.size() ; j++) {
				center[j] = floor(center[j]/gridSize);
				cen_tmp[j] = floor(cen_tmp[j]/gridSize);
			}
			int flag = 0;
			for (int j = 0 ; j < center.size() ; j++) {
				flag += abs(cen_tmp[j] - center[j]);
			}
			if (flag == 1 || flag == 0) {
				neighbour.push_back(i);
			}
		}

		return neighbour;
	}

	private:void sortNeighbour (vector<int>& neighbour,Clustering<MiniCluster> x,MiniCluster mc) {
		vector<double> distance;
		for (int i = 0 ; i < neighbour.size() ; i++) {
			distance.push_back(dist(x.get(neighbour[i]).getcenter(), mc.getcenter()));
		}
		quick_sort(distance, 0, distance.size()-1, neighbour);
	}
	
	private:void expandCluster (Clustering<MiniCluster> *cluster,Clustering<MiniCluster> x,MiniCluster core,vector<int> neighbour,vector<double> stat) {
		int size = neighbour.size();
		for (int i = 0 ; i < size ; i++) {
			if (!cl[neighbour[i]]) {
				cl[neighbour[i]] = 1;
				MiniCluster mc = x.get(neighbour[i]);
				cluster->add(mc);
				vector<int> neighbour2 = cmcNeighbourhood(x,mc);
				if (neighbour2.size() >= minpts) {
					vector<double> stat2 = coreNeighbourStat(neighbour2, x, mc);
					sortNeighbour(neighbour2, x, mc);
					MiniCluster nearestNeighbour = x.get(neighbour2[0]);
					vector<int> neighbour3 = cmcNeighbourhood(x, nearestNeighbour);
					if (neighbour3.size() >= minpts) {
						if (stat2[0] >= (stat[0] - stat[1]) && stat2[0] <= (stat[0] + stat[1])) {
							for (int j = 0 ; j < neighbour3.size() ; j++) {
								if(find(neighbour.begin(), neighbour.end(), neighbour3[j]) == neighbour.end()){
									neighbour.push_back(neighbour3[j]);
								}
							}
							size = neighbour.size();
							stat = coreNeighbourStat(neighbour, x, core);
						}
					}
				}
			}
		}	
	}

	private:vector<double> coreNeighbourStat (vector<int> neighbour,Clustering<MiniCluster> x,MiniCluster cmcp) {
		vector<double> stat;
		stat.push_back(0.0);
		stat.push_back(0.0);

		for (int i = 0 ; i < neighbour.size() ; i++) {
			MiniCluster mc = x.get(neighbour[i]);
			stat[0] += dist(cmcp.getcenter(), mc.getcenter());
		}
		stat[0] /= neighbour.size();
		for (int i = 0 ; i < neighbour.size() ; i++) {
			MiniCluster mc = x.get(neighbour[i]);
			stat[1] += pow(dist(cmcp.getcenter(), mc.getcenter()) - stat[0],2);
		}
		stat[1] /= neighbour.size();

		return stat;
	}

	public:List result() {
		int N = 0;
		cout << "\nForming Final Cluster !!!\n";
		cl = new int[cmc.sizes()];
		vector<Clustering<MiniCluster> > cluster;	
		for (int i = 0 ; i < cmc.sizes() ; i++) {
			cl[i] = 0;
		}
		for (int i = 0 ; i < cmc.sizes() ; i++) {
			if (cl[i] == 0) {
				cl[i] = 1;
				MiniCluster mc = cmc.get(i);
				vector<int> neighbour = cmcNeighbourhood(cmc,mc);
				if (neighbour.size() >= minpts) {
					cluster.push_back(Clustering<MiniCluster>());
					cluster[cluster.size() - 1].add(mc);
					vector<double> stat = coreNeighbourStat(neighbour, cmc, mc);
					sortNeighbour(neighbour, cmc, mc);
					MiniCluster nearestNeighbour = cmc.get(neighbour[0]);
					vector<int> neighbour2 = cmcNeighbourhood(cmc,nearestNeighbour);
					if (neighbour2.size() >= minpts) {
						expandCluster(&cluster[cluster.size() - 1], cmc, mc, neighbour, stat);
						N += cluster[cluster.size() - 1].sizes();
					}
				} else {
				}	
			}	
		}
		cout << "Number of micro-cluster: " << cmc.sizes() << endl;
		cout << "Number of marcro-cluster: "<< cluster.size() << endl;
		// cout << "gridSize = " << grid_list.size() << endl;
		List centroid(N);
		int n = 0;
		for (int i = 0 ; i < cluster.size() ; i++) {
			for (int j = 0 ; j < cluster[i].sizes() ; j++) {
				Clustering<MiniCluster> clus = cluster[i];
				std::vector<double> v = clus.get(j).getcenter();
				v.push_back(double(i));
				centroid[n++] = v;
			}
		}
		return centroid;
	}
};


RCPP_MODULE(MOD_Mudi){
  
  class_<Mudi>("Mudi")
  .constructor<int, double, double, double, int>()
  .method("trainOnInstance", &Mudi::trainOnInstance)
  .method("result", &Mudi::result)
  ;
}