#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

vector<double> composition(vector < vector < double>> A, vector<double> b) {
	vector<double> comp(17, 0);
	for (int i = 0; i < 17; i++)
		for (int j = 0; j < 17; j++) {
			comp[i] += A[i][j] * b[j];
		}
	return comp;
}

vector<double> subtraction(vector<double> x, vector<double> y) {
	vector<double> sub(17);
	for (int i = 0; i < 17; i++) {
		sub[i] = x[i] - y[i];
	}
	return sub;
}

double composition(vector<double> x, vector<double> y) {
	double ans = 0;
	for (int i = 0; i < x.size(); i++) {
		ans += x[i] * y[i];
	}
	return ans;
}
vector<double> composition(vector<double> x, double y) {
	vector<double> comp = x;
	for (int i = 0; i < x.size(); i++) {
		comp[i] *= y;
	}
	return comp;
}

vector<double> addition(vector<double> x, vector<double>y) {
	vector<double> add = x;
	for (int i = 0; i < x.size(); i++) {
		add[i] += y[i];
	}
	return add;
}

double norm(vector<double> A)
{
	double sum = 0;
	for (auto u : A) sum += u * u;
	return sqrt(sum);
}

void my_pcg(vector < vector < double>>& A, vector<double>& b, double precision = false, int iter_limit = false) {
	int iters = 0;

	vector<double> x = b;
	vector<double> r = subtraction(b, composition(A, b));
	vector<double> z = r;

	bool conditions_met = true;

	while (conditions_met) {
		iters++;

		if (composition(composition(A, z), z) <= 0 || isinf(composition(composition(A, z), z))) {
			cout << "Method discrepancy.\n";
			return;
		}

		double alpha = composition(r, r) / composition(composition(A, z), z);

		x = addition(x, composition(z, alpha));
		vector<double> r_prev = r;
		r = subtraction(r, composition(composition(A, z), alpha));

		if (composition(r_prev, r_prev) == 0 || isinf(composition(r_prev, r_prev))) {
			cout << "Method discrepancy.\n";
			return;
		}

		double beta = composition(r, r) / composition(r_prev, r_prev);

		z = addition(r, composition(z, beta));

		if (iter_limit) {
			conditions_met *= iters < iter_limit;
		}
		if (precision) {
			conditions_met *= norm(r) > precision;
		}
	}
	if (iters >= iter_limit) {
		cout << "There is no accuracy.\n";
	}
	cout << iters << '\n';
	cout << "x = \n";
	for (auto u : x) {
		cout << u << ' ';
	}
	cout << '\n';
}

int main() {
	cout << "Ax = b\n\n";
	vector<vector<double> > A;
	A.resize(17, vector<double>(17));
	//Ax=b
	cout << "A = \n";
	for (int i = 0; i < 17; i++) {
		for (int j = 0; j < 17; j++) {
			cin >> A[i][j];
		}
	}
	cout << '\n';

	vector<double> b(17);
	cout << "b = \n";
	for (int i = 0; i < 17; i++) {
		cin >> b[i];
	}
	cout << "\n";

	my_pcg(A, b, 0.000001, 50);
	return 0;
}