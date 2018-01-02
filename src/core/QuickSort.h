#include <Rcpp.h>
using namespace Rcpp;
using namespace std;


void quick_sort (vector<double>& arr, int left, int right,vector<int>& relation) {
	int i = left;
	int j = right;
	double pivot = arr[(left + right)/2];
	double tmp;
	int tmp2;
	while (i <= j) {
		while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;

        if (i <= j) {
            tmp = arr[i];
            arr[i] = arr[j];
            arr[j] = tmp;

            /*Swap relation*/
            tmp2 = relation[i];
            relation[i] = relation[j];
            relation[j] = tmp2;

            i++;
            j--;

            
        }
	}
	if (left < j)
       quick_sort(arr, left, j, relation);
    if (i < right)
       quick_sort(arr, i, right, relation);

}