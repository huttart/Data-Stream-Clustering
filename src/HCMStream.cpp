#include"core/SphereCluster.h"
#include"core/Clustering.h"
#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

double dotProduct (vector<double> x1,vector<double> x2) {
	double sum = 0;
	for (int i = 0 ; i < x1.size() ; i++) {
		sum += x1[i]*x2[i];
	}
	return sum;
}

double vectorSize (vector<double> x) {
	double size = 0;
	for (int i = 0 ; i < x.size() ; i++) {
		size += pow(x[i],2);
	}
	return sqrt(size);
}

vector<double> unitVector (vector<double> x) {
	double size = vectorSize(x);
	for (int i = 0 ; i < x.size() ; i++) {
		x[i] /= size;
	}
	return x;
}

vector<double> plusVector (vector<double> x1,vector<double> x2) {
	vector<double> plus;
	for (int i = 0 ; i < x1.size() ; i++) {
		plus.push_back(x1[i]+x2[i]);
	}
	return plus;
}

vector<double> minusVector (vector<double> x1,vector<double> x2) {
	vector<double> plus;
	for (int i = 0 ; i < x1.size() ; i++) {
		plus.push_back(x1[i]-x2[i]);
	}
	return plus;
}

vector<double> scalarProduct (vector<double> x1,double scalar) {
	vector<double> plus;
	for (int i = 0 ; i < x1.size() ; i++) {
		plus.push_back(x1[i]*scalar);
	}
	return plus;
}

vector<double> scalarDivide (vector<double> x1,double scalar) {
	vector<double> plus;
	for (int i = 0 ; i < x1.size() ; i++) {
		plus.push_back(x1[i]/scalar);
	}
	return plus;
}



class SphereMicro: public Sphere {
public:
	SphereMicro (double r,vector<double> point) : Sphere(r,point) {
		// this->r = r;
		// this->center = point;
		// N = 1;
	}

	public:void insert (vector<double> point) {
		this->center = scalarDivide(plusVector(scalarProduct(this->center ,this->N), point), this->N+1);
		this->N++;
	}
};

class Denpoint {
public:
	bool covered;

	Denpoint () {
		this->covered = false;
	}

	public:void set_covered () {
		this->covered = true;
	}
};

class CylinderMicro {
public:
	int N;

	double L;
	double r;

	vector<double> I;
	vector<double> center;

	CylinderMicro (double r,vector<SphereMicro> x) {
		vector<double> center_tmp;
		vector<vector<double> > C;
		this->N = 0;
		this->r = r;
		for (int i = 0 ; i < x.size() ; i++) {
			this->N++;
			if (i == 0) {
				center_tmp = x[i].center;
			} else {
				for (int j = 0 ; j < center_tmp.size() ; j++) {
					center_tmp[j] += x[i].center[j];
				}
			}
		}
		for (int j = 0 ; j < center_tmp.size() ; j++) {
			this->center.push_back(center_tmp[j]/x.size());
		}
		for (int i = 0 ; i < x.size() ; i++) {
			vector<double> tmp;
			for (int j = 0 ; j < this->center.size() ; j++) {
				tmp.push_back(x[i].center[j] - this->center[j]);
			}
			C.push_back(tmp);
		}
		if(x.size() == 2) {
			vector<double> c1 = x[0].center;
			vector<double> c2 = x[1].center;
			double vector_size = 0;
			for (int j = 0 ; j < c1.size() ; j++) {
				I.push_back(c1[j] - c2[j]);
				vector_size += pow(c1[j] - c2[j],2);
			}
			vector_size = sqrt(vector_size);
			for (int j = 0 ; j < I.size() ; j++) {
				I[j] /= vector_size;
			}
		} else {
			double maxdis = 10000000;
			int index = 0;
			for (int i = 0 ; i < C.size() ; i++) {
				double dis = 0;
				for (int j = 0 ; j < C[i].size() ; j++) {
					dis += pow(C[i][j],2);
				}
				dis = sqrt(dis);
				if (dis > maxdis) {
					maxdis = dis;
					index = i;
				}
			}
			double vector_size = 0;
			for (int i = 0 ; i < C[index].size() ; i++) {
				I.push_back(C[index][i]);
				vector_size += pow(C[index][i],2);
			}
			vector_size = sqrt(vector_size);
			for (int j = 0 ; j < I.size() ; j++) {
				I[j] /= vector_size;
			}
		}

		double maxdis = 10000000;
		int index = 0;
		for (int i = 0 ; i < C.size() ; i++) {
			double dis = 0;
			for (int j = 0 ; j < C[i].size() ; j++) {
				dis += pow(C[i][j],2);
			}
			dis = sqrt(dis);
			if (dis > maxdis) {
				maxdis = dis;
				index = i;
			}
		}
		for (int i = 0 ; i < C[index].size() ; i++) {
			L += C[index][i] * I[i];
		}

		L += r;

	}

	public:bool mergeWithSM (SphereMicro *x,double dpara,double length) {
		this->N += x->N;
		if (dpara+r > this->L) {
			this->L = (dpara+r+this->L)/2;
			if(this->L > length) {return false;}
			double tmp = (dotProduct(this->I,minusVector(x->center,this->center)));
			vector<double> tmp2 = scalarProduct(this->I,tmp);
			tmp2 = scalarProduct(tmp2,(dpara-this->L+r)/2);
			this->center = plusVector(this->center,tmp2);
		}
		return true;
	}

};

class HCMStream {
private:
	double r;
	int np;
	int numberOfCluster;
	double length;

	Clustering<SphereMicro> SM;
	Clustering<SphereMicro> SMTemp;
	Clustering<CylinderMicro> CM;

	vector<vector<int> > list;
	// vector<vector<int> > cl;
	int *cl;

public:
	HCMStream (double r,int np,double length) {
		this->r = r;
		this->np = np;
		this->length = length;
		this->numberOfCluster = 0;
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

	private:int nearestCM(vector<double> input,Clustering<CylinderMicro> x) {
		int min = -1;
		double distance = 1000000;
		int size = x.sizes();
		for (int i = 0 ; i < size ; i++) {
			CylinderMicro mc = x.get(i);
			double dpara = (dotProduct(minusVector(input,mc.center), mc.I));
			vector<double> tmp = scalarProduct(mc.I,dpara);
			dpara = abs(dpara);
			tmp = minusVector(minusVector(input, mc.center), tmp);
			double dperp = vectorSize(tmp);
			// cout << "dpara = " << dpara << endl;
			// cout << "dperp = " << dperp << endl;
			if (dpara <= mc.L && dperp <= r) {
				if (dperp < distance) {
					distance = dperp;
					min = i;
				}
			}
		}
		return min;	
	}

	private:int nearestSM(vector<double> input,Clustering<SphereMicro> x) {
		int min = -1;
		double distance = 1000000;
		int size = x.sizes();
		for (int i = 0 ; i < size ; i++) {
			SphereMicro mc = x.get(i);
			double dis = vectorSize(minusVector(input,mc.center));
			if (dis <= r) {
				if (dis < distance) {
					distance = dis;
					min = i;
				}
			}
		}
		return min;	
	}

	private:int mergeSMtoCM (SphereMicro *y,Clustering<CylinderMicro> *x) {
		int min = -1;
		bool merge = false;
		double distance = 1000000;
		int size = (*x).sizes();
		double dpararell;
		vector<double> input = y->center;
		for (int i = 0 ; i < size ; i++) {
			CylinderMicro mc = (*x).get(i);
			double dpara = (dotProduct(minusVector(input,mc.center), mc.I));
			vector<double> tmp = scalarProduct(mc.I,dpara);
			dpara = abs(dpara);
			tmp = minusVector(minusVector(input, mc.center), tmp);
			double dperp = vectorSize(tmp);
			// cout << "dpara = " << dpara << endl;
			// cout << "dperp = " << dperp << endl;
			if (dpara <= mc.L+r && dperp <= r*2) {
				if (dperp < distance) {
					distance = dperp;
					min = i;
					dpararell = dpara;
				}
			}
		}
		if (min >= 0) {
			CylinderMicro *mc = (*x).getreal(min);	
			if(!mc->mergeWithSM(y,dpararell,this->length)){
				min = -1;
			}
		}
		return min;
	}

	private:bool creatNewCylinder (SphereMicro mc,Clustering<SphereMicro> *x) {
		vector<int> sMerge_index;
		vector<SphereMicro> sMerge;
		for (int i = 0 ; i < x->sizes() ; i++) {
			SphereMicro x2 = x->get(i);
			if (dist(x2.center,mc.center) <= 2*r && x2.N >= np) {
				sMerge_index.push_back(i);
				sMerge.push_back(x2);
			}
		}
		if (sMerge_index.size() >= 2) {

			CylinderMicro newCylinder(r,sMerge);
			CM.add(newCylinder);
			for (int i = sMerge_index.size()-1 ; i >= 0 ; i--) {
				x->remove(sMerge_index[i]);
			}
			return true;
		}

		return false;


	}

	public:void trainOnInstance(NumericVector point) {
		vector<double> input = as<vector<double> >(point);
		bool merge = false;
		int nearest;
		if(CM.sizes() > 0) {
			nearest = nearestCM(input,CM);
			if (nearest >= 0) {
				CM.getreal(nearest)->N++;
				merge = true;
			}
		}

		if (!merge) {
			nearest = nearestSM(input,SM);
			if (nearest >= 0) {
				SphereMicro *mc = SM.getreal(nearest);
				mc->insert(input);
				merge = true;
				if (mc->N > this->np) {
					int merge2 = mergeSMtoCM(mc,&CM);
					if (merge2 >= 0) {
						SM.remove(nearest);
					} else {
						creatNewCylinder(SM.get(nearest),&SM);
					}

				}
			}
		}

		if (!merge) {
			SM.add(SphereMicro(r,input));
		}
	}

	private:vector<int> getNeighbourhoodIDs(vector<double> input,Clustering<SphereMicro> x) {
		vector<int> neighbourhood;
		int size = x.sizes();
		for (int i = 0 ; i < size ; i++) {
			SphereMicro mc = x.get(i);
			double dis = vectorSize(minusVector(input,mc.center));
			cout << "dis = " << dis << endl;

			if (dis <= 2*r && dis != 0) {
				neighbourhood.push_back(i);
			}
		}
		return neighbourhood;	
	}

	// private:void expandCluster (vector<Denpoint> *cover,Clustering<SphereMicro> sm,SphereMicro mc) {
	// 	for (int i = 0 ; i < sm.sizes() ; i++) {
	// 		SphereMicro mc2 = sm.get(i);
	// 		if (dist(mc2.center,mc.center) < 2*this->r && !(*cover)[i].covered) {
	// 			cout << "i = " << i << endl;
	// 			cout << "covered = " << (*cover)[i].covered << endl;
	// 			vector<int> tmp;
	// 			tmp.push_back(numberOfCluster);
	// 			tmp.push_back(list[i][1]);
	// 			cl.push_back(tmp);
	// 			(*cover)[i].set_covered();
	// 			expandCluster(cover,sm,mc2);
	// 		}
	// 	}
	// }

	public:List finalCluster(const vector<double> center = vector<double>(),double eps = 0,bool starting = false) {
		if (starting && SM.sizes() > 0) {
			cl = new int[SM.sizes()];
			for (int i = 0 ; i < SM.sizes() ; i++) {
				cl[i] = 0;
			}
			for (int i = 0 ; i < SM.sizes() ; i++) {
				SphereMicro mc = SM.get(i);
				if (!cl[i]) {
					numberOfCluster += 1;
					cl[i] = numberOfCluster;
					finalCluster(mc.center,eps);
				}
			}
		} else {
			for (int i = 0 ; i < SM.sizes() ; i++) {
				SphereMicro mc = SM.get(i);
				if (!cl[i]) {
					vector<double> center2 = mc.center;
					double distance = dist(center,center2);
					if (distance <= eps*2) {
						cl[i] = numberOfCluster;
						finalCluster(center2,eps);
					}
				}
			}
		}
	}

	private:void makeFinalCluster () {
		for (int i = 0 ; i < SM.sizes() ; i++) {
			SphereMicro sm = SM.get(i);
			if (sm.N <= this->np) {
				SM.remove(i);
			} 
		}
		for (int i = 0 ; i < CM.sizes() ; i++) {
			CylinderMicro cm = CM.get(i);
			SM.add(SphereMicro(r,cm.center));
			SM.getreal(SM.sizes()-1)->N = 1000000;
			vector<int> tmp;
			tmp.push_back(SM.sizes()-1);
			tmp.push_back(i);
			list.push_back(tmp);
			vector<double> I = cm.I;
			vector<double> currentCenter = cm.center;
			vector<double> newCenter = plusVector(scalarProduct(I,r),currentCenter);
			while (true) {
				SM.add(SphereMicro(r,newCenter));
				SM.getreal(SM.sizes()-1)->N = 100000;
				vector<int> tmp2;
				tmp2.push_back(SM.sizes()-1);
				tmp2.push_back(i);
				list.push_back(tmp2);
				currentCenter = newCenter;
				if (dist(currentCenter,cm.center)+r >= cm.L) {
					break;
				}
				newCenter = plusVector(scalarProduct(I,r),currentCenter);
			}
			I = scalarProduct(I,-1);
			currentCenter = cm.center;
			newCenter = plusVector(scalarProduct(I,r),currentCenter);
			while (true) {
				SM.add(SphereMicro(r,newCenter));
				SM.getreal(SM.sizes()-1)->N = 100000;
				vector<int> tmp3;
				tmp3.push_back(SM.sizes()-1);
				tmp3.push_back(i);
				list.push_back(tmp3);
				currentCenter = newCenter;
				if (dist(currentCenter,cm.center)+r >= cm.L) {
					break;
				}
				newCenter = plusVector(scalarProduct(I,r),currentCenter);
			}
		}
		finalCluster(vector<double>(),this->r,true);

	}

	public:List result () {
		makeFinalCluster();
		int cm_no = 0;
		vector<int> clusterOfCm;
		for (int i = 0 ; i < list.size() ; i++) {
			if (list[i][1] == cm_no) {
				clusterOfCm.push_back(cl[list[i][0]]);
				cm_no++;
			}
		}
		List x(CM.sizes()+SM.sizes());
		int n = 0;
		for (int i = 0 ; i < CM.sizes() ; i++) {
			CylinderMicro cm = CM.get(i);
			vector<double> tmp = cm.center;
			for (int j = 0 ; j < cm.I.size() ; j++) {
				tmp.push_back(cm.I[j]);
			}
			tmp.push_back(cm.L);
			tmp.push_back(cm.N);
			tmp.push_back(clusterOfCm[i]);
			x[n++] = tmp;
		}

		for (int i = 0 ; i < SM.sizes() ; i++) {
			vector<double> tmp = SM.get(i).center;
			tmp.push_back(SM.get(i).N);
			x[n++] = tmp;
		}

		// for (int i = 0 ; i < SMTemp.sizes() ; i++) {
		// 	vector<double> tmp = SMTemp.get(i).center;
		// 	tmp.push_back(2.0);
		// 	x[n++] = tmp;
		// }

		return x;

	}
};


RCPP_MODULE(MOD_HCM){
  
  class_<HCMStream>("HCMStream")
  .constructor<double, int, double>()
  .method("trainOnInstance", &HCMStream::trainOnInstance)
  .method("result", &HCMStream::result)
  ;
}