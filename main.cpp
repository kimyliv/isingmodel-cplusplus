#include <iostream>
#include <random>
#include <vector>

namespace MCnameSpace {
    class MonteCarloSimulation {

    private:
        int L; // lattice size
        int J; // couplingconst
        double Tmin;
        double Tmax;
        double stepsize;
        std::vector<std::vector<int>> lattice;
        std::random_device rd;
        std::mt19937 gen;
        std::uniform_real_distribution<> dis;
        std::uniform_int_distribution<> disint;
        std::uniform_int_distribution<> initpos;


    public:
        MonteCarloSimulation(int size, int couplingconst, double Tempmin, double Tempmax, double Tempstep) : L(size),
                               J(couplingconst), Tmin(Tempmin), Tmax(Tempmax),stepsize(Tempstep),gen(rd()),dis(0.0,1.0),
                               disint(-1.0,1.0),initpos(0.0,L -1) {
        };

        void Initialcond() {
            lattice.resize(L, std::vector<int>(L, 1));
        }

        void thermalization(double t){
            int x = initpos(gen);
            int y = initpos(gen);
            double thermalizationsteps = 50*L*L; // 50 mcsteps for thermalization

            for (int i = 0; i < thermalizationsteps; ++i) {

                x = x + disint(gen);
                y = y + disint(gen);

                double dE = deltaE(&x, &y);
                //std::cout<<"Random number: " << dis(gen) << " exp value: " << exp(-dE/t) << std::endl;

                if (dE < 0 || dis(gen) < exp(-dE / t)) {
                    lattice[x][y] *= -1;
                }
            }
        }

        double deltaE(int *ip, int *jp) {

            int len = L - 1;
            int i = *ip;
            int j = *jp;

            if (i > len) {
                i = 0;
            }
            if (j > len) {
                j = 0;
            }
            if (i < 0) {
                i = len;
            }
            if (j < 0) {
                j = len;
            }

            int left = (i - 1 + len) % len;
            int right = (i + 1) % len;
            int up = (j - 1 + len) % len;
            int down = (j + 1) % len;

            double dE = 2 * J * (lattice[i][j]) *
                        (lattice[left][j] + lattice[right][j] + lattice[i][down] + lattice[i][up]);
            *jp = j;
            *ip = i;
            return dE;

        }

        void Mcstep(double t) {

            int x = initpos(gen);
            int y = initpos(gen);

            for (int i = 0; i < L * L; ++i) {

                x = x + disint(gen);
                y = y + disint(gen);

                double dE = deltaE(&x, &y);

                if (dE < 0 || dis(gen) < exp(-dE / t)) {
                    lattice[x][y] *= -1;
                }
            }
        }

        double magnetization() {

            double m = 0;

            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {
                    m += lattice[i][j];
                }
            }
            return abs(m) / (L * L);
        }

        double Energycalculations(){
            int len = L-1;
            double energy = 0.0;
            for (int i = 0; i < L; ++i) {
                for (int j = 0; j < L; ++j) {

                    int left = (i - 1 + len) % len;
                    int right = (i + 1) % len;
                    int up = (j - 1 + len) % len;
                    int down = (j + 1) % len;

                    energy -= J*lattice[i][j]*(lattice[left][j] + lattice[right][j] + lattice[i][up] + lattice[i][down]);

                }

            }
            return energy /(L*L);
        }

        void simulation(int McSteps) {
            double Tempval,Tempvale;
            double M, M2, M4,gbinder, suscept, specheat;
            double En, En2;

            std::cout<<" -Temperature- "<<" -Magnetization- "<< " -Energy- "<< " -Specific heat- "
            << "-Magnetic susceptibility- " << " -Binder cumulant- "<<std::endl;
            for (double t = Tmin; t < Tmax; t = t + stepsize) {

                M=M2=M4=0;
                En = En2 = 0;
                Initialcond();
                thermalization(t);

                for (int i = 0; i < McSteps; ++i) {

                    Mcstep(t);
                    Tempval = magnetization();
                    Tempvale = Energycalculations();
                    En += Tempvale;
                    En2 += Tempvale*Tempvale;

                    M += Tempval;
                    M2 += Tempval*Tempval;
                    M4 += Tempval*Tempval*Tempval*Tempval;
                }

                En = En/(McSteps);
                En2 = En2/(McSteps);
                M = M/(McSteps);
                M2 = M2/(McSteps);
                M4 = M4/(McSteps);

                specheat = (En2 - En*En)/(L*L*t*t);
                suscept = (L*L)*(M2 - M*M)/t;
                gbinder = 1- M4/(3*M2*M2);

                std::cout <<"    "<< t <<"          "<< M <<"          "<< En <<"          "
                << specheat <<"          "<< suscept <<"          "<< gbinder <<"          "<< std::endl;
            }
        }

    };
}
int main() {

    double tmin = 0.5;
    double tmax = 5.0;
    double tstep = 0.1;

    int gridsize = 20;
    int J = 1;

    MCnameSpace::MonteCarloSimulation model(gridsize,J,tmin,tmax,tstep);
    model.simulation(10000);

    return 0;
}
