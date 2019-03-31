#include<cmath>
#include<vector>
#include<ctime>
#include<random>
#include<iostream>
#define pi 3.1415926535
#define E  2.71828182845904523536

using namespace std;

double randDouble(double min, double max)
{
	static default_random_engine engine(time(nullptr));
	uniform_real_distribution<double> dis(min, max);
	return dis(engine);
}
class Particle
{
public:
	vector<double> position;
	vector<double>velocity;
	vector<double>pBest;
	double Pc;
	vector<int> fi;
	int times;
};

bool better(double a, double b)
{
	if (a < b)
		return true;
	else
		return false;
}

class PSO
{
public:

	PSO(int dim, int m, int Tmax, double max, double min, double c, 
		double wmax, double wmin, double dt, double percent,int R)
	{
		this->R = R;
		this->dim = dim;
		this->m = m;
		this->Tmax = Tmax;
		this->max = max;
		this->min = min;
		this->c = c;
		this->wmax = wmax;
		this->wmin = wmin;
		this->dt = dt;
		this->percent = percent;
		particles.resize(m);
	}


	bool inRange(vector<double> pos)
	{
		for (int i = 0; i < pos.size(); i++)
		{
			if (pos[i] > max)
				return false;
			if (pos[i] < min)
				return false;
		}
		return true;
	}

	double fitnessFunction(vector<double>pos)
	{
		double result = 0.0;

		for (int i = 0; i < dim; i++)
		{
			double temp = 0.0;
			for (int j = 0; j < i; j++)
			{
				temp += pos[j];
			}
			result += pow(temp, 2);
		}
		return result;
	}

	void initialParticles(int i)
	{
		particles[i].position.resize(dim);
		particles[i].velocity.resize(dim);
		particles[i].pBest.resize(dim);
		particles[i].fi.resize(dim);
		particles[i].Pc = 0.05 + 0.45*(exp(10 * (i - 1) / (m - 1)) - 1) / (exp(10) - 1);
		particles[i].times = 0;
		for (int j = 0; j < dim; j++)
		{
			double range = percent * (max - min);
			particles[i].position[j] = randDouble(this->min, this->max);
			particles[i].velocity[j] = randDouble(-range, range);
			particles[i].pBest[j] = particles[i].position[j];
			particles[i].fi[j] = i;
		}
		
	}

	void initialAllParticles()
	{
		initialParticles(0);

		for (int i = 1; i < m; i++)
		{
			initialParticles(i);
		}
	}

	void inertiaWeight()
	{
		//w = randDouble(0.4, 0.6);
		double t = T / ((double)Tmax);
		w = wmax - (wmax - wmin)*t;
	}

	void updateParticle(int i)
	{
		static default_random_engine engine(time(nullptr));
		uniform_real_distribution<double> dis(0, m - 1);

		if(particles[i].times>=R)
		{
			particles[i].times = 0;
			for (int j = 0; j < dim; j++)
			{
				particles[i].fi[j] = i;
				if (randDouble(0, 1) < particles[i].Pc)
				{
					
					int a = dis(engine);
					int b = dis(engine);
					while (1)
					{
						a = dis(engine);
						b = dis(engine);
						if (a != b && a != i && b != i)
							break;
					}
					if (fitnessFunction(particles[b].pBest) < fitnessFunction(particles[a].pBest))
						particles[i].fi[j] = b;
					else
						particles[i].fi[j] = a;
					
					cout << a << " " << b << endl;
				}
				
			}
			
		}
		

		for (int j = 0; j < dim; j++)
		{
			double range = percent * (max - min);

			particles[i].velocity[j] = w * particles[i].velocity[j] +
				c * randDouble(0, 1) * (particles[particles[i].fi[j]].pBest[j] - particles[i].position[j]);
			if (particles[i].velocity[j] > range)
				particles[i].velocity[j] = range;

			if (particles[i].velocity[j] < -range)
				particles[i].velocity[j] = -range;
			particles[i].position[j] += dt * particles[i].velocity[j];
		}

		

		if (inRange(particles[i].position))
		{
			if (fitnessFunction(particles[i].position) < fitnessFunction(particles[i].pBest))
			{
				particles[i].pBest = particles[i].position;
				particles[i].times = 0;
			}
			else
				particles[i].times++;
		}

	}


	void updateAllParticles()
	{
		inertiaWeight();
		for (int i = 0; i < m; i++)
		{
			updateParticle(i);
		}
		T++;
	}

	double getFitness()
	{
		int index = 0;
		for (int i = 0; i < m; i++)
		{
			if (fitnessFunction(particles[i].pBest) < fitnessFunction(particles[index].pBest))
				index = i;
		}
		return fitnessFunction(particles[index].pBest);
	}
private:
	int dim;
	int m;//number of instances

	int T;
	int Tmax;

	double w;
	double max;
	double min;
	double c;
	double wmax;
	double wmin;
	int R;

	double dt;//时间步长
	double percent;


	vector<Particle> particles;


};

void run(vector<double>& result1)
{
	int dim = 30;
	int m = 20;
	int Tmax = 2000;
	double max = 100;
	double min = -100;
	double c = 1.49445;
	double wmax = 0.9;
	double wmin = 0.4;
	double dt = 1.0;
	double percent = 0.2;
	int R = 7;

	PSO pso = PSO(dim, m, Tmax, max, min, c, wmax, wmin, dt, percent, R);
	pso.initialAllParticles();

	vector<double>fitness;
	fitness.push_back(pso.getFitness());

	for (int i = 0; i < Tmax; i++)
	{
		pso.updateAllParticles();
		cout << ":";
		//fitness.push_back(pso.getFitness());
		fitness.push_back(pso.getFitness());
		cout << "第" << i << "次迭代结果：";
		cout << ", fitness = " << pso.getFitness() << endl;



	}
	result1 = fitness;
}

int main()
{

	int times = 5;
	int interval = 10;
	vector<double> result1;

	run(result1);

	for (int i = 1; i < times; i++)
	{
		vector<double> result1_temp;
		run(result1_temp);
		for (int j = 0; j < result1_temp.size(); j++)
		{
			result1[j] += result1_temp[j];
		}
	}
	for (int j = 0; j < result1.size(); j++)
	{
		result1[j] /= times;
	}

	for (int j = 0; j < result1.size(); j++)
	{
		if (j%interval == 0)
			cout << result1[j] << " ";
	}

	system("pause");
}