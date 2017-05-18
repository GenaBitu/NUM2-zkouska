using namespace std;

#include <iostream>
#include <fstream>
#include <cmath>
#include <tuple>
#include <functional>

const float PI = 3.14159265352;

void write_solution(float c_1_exact, float c_2_exact, float c_1_approx, float c_2_approx)
{
	string head{R"(set terminal qt size 1000,1000 enhanced font 'Verdana,10' persist
# Line width of the axes
set border linewidth 1.5
# Line styles
set style line 1 linecolor rgb '#0060ad' linetype 1 linewidth 2
set style line 2 linecolor rgb '#dd181f' linetype 1 linewidth 2
# Axes label
set xlabel 'x'
set ylabel 'y'
# Axes ranges
set xrange [0:0.5*pi]
set xtics ('0' 0, 'π/8' 0.125 * pi, 'π/4' 0.25 * pi, '3π/8' 0.375 * pi, 'π/2' 0.5 * pi)
# Paramters
c_1_exact = )"};

    string line1{R"(
c_2_exact = )"};

    string line2{R"(
c_1_approx = )"};

    string line3{R"(
c_2_approx = )"};

    string tail{R"(
# Fuction
exact(x) = c_1_exact * exp(4 * x) + c_2_exact * exp(-4 * x) + 0.5
approx(x) = c_1_approx * exp(4 * x) + c_2_approx * exp(-4 * x) + 0.5
# Plot
plot exact(x) title 'exaktní řešení' with lines linestyle 1, \
     approx(x) title 'přibližné řešení' with lines linestyle 2)"};

	ofstream file{"./solution.gnu"};

	file << head << c_1_exact << line1 << c_2_exact << line2 << c_1_approx << line3 << c_2_approx << tail;
}

tuple<float, float> halve(float alpha, float beta, float epsilon)
{
	function<float(float, float)> c_1{[](float a, float c)->float{return 0.5 * a + 0.125 * c - 0.25;}};
	function<float(float, float)> c_2{[](float a, float c)->float{return 0.5 * a - 0.125 * c + 0.25;}};

	function<float(float)> F{[alpha, beta, c_1, c_2](float c)->float{return c_1(alpha, c) * exp(2.0 * PI) + c_2(alpha, c) * exp(-2.0 * PI) + 0.5 - beta;}};

	float gamma_left, gamma_right;

	for(int i = 0; i < 100000; ++i)
	{
		gamma_left = static_cast<float>(rand()) / static_cast<float>(RAND_MAX/2000) - 1000;
		gamma_right = static_cast<float>(rand()) / static_cast<float>(RAND_MAX/2000) - 1000;
		if(F(gamma_left) * F(gamma_right) < 0.0)
		{
			break;
		}
	}
	if(F(gamma_left) * F(gamma_right) >= 0.0)
	{
		cout << "Nepodařilo se najít počáteční stav metody půlení intervalů" << endl;
		abort();
	}

	cout << "Půlím..." << endl;

	while(true)
	{
		float gamma = 0.5 * (gamma_left + gamma_right);
		float middle = F(gamma);
		if(abs(middle) < epsilon)
		{
			return make_tuple(c_1(alpha, gamma), c_2(alpha, gamma));
		}
		if(middle * F(gamma_right) > 0.0)
		{
			gamma_right = gamma;
		}
		else
		{
			gamma_left = gamma;
		}
	}
}

int main()
{
	float alpha;
	float beta;
	float epsilon;
	cout << "Zadejte parametr alfa: " << endl;
	cin >> alpha;
	cout << "Zadejte parametr beta: " << endl;
	cin >> beta;
	cout << "Zadejte parametr epsilon: " << endl;
	cin >> epsilon;

	float c_1_exact{((beta + 0.5f) + exp(-2.0f * PI) * (0.5f - alpha)) / (2.0f * sinh(2.0f * PI))};
	float c_2_exact{alpha - 0.5f - c_1_exact};

	tuple<float, float> c_approx = halve(alpha, beta, epsilon);

	write_solution(c_1_exact, c_2_exact, get<0>(c_approx), get<1>(c_approx));
}