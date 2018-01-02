#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

template <class T>
class Clustering {
private:
	vector<T> cluster;
	int size;
public:
	Clustering () {
		this->size = 0;
	}

	public:void add(T x) {
		cluster.push_back(x);
		this->size += 1;
	}

	public:void remove(int index) {
		cluster.erase(cluster.begin() + index);
		this->size -= 1;
	}

	public:T get(int index) {
		return cluster[index];
	}

	public:T* getreal(int index) {
		return &cluster[index];
	}

	public:int sizes() {
		return this->size;
	}
};



