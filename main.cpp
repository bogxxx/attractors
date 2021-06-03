#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <vector>
#include <omp.h>
#include <ctime>
#include <string>
#include <cmath>
#include <sstream>
#include <cstdlib>

using namespace std;
using namespace boost::numeric::odeint;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

const double xb = -5.0;
const double xe = 5.0;
const double dx = 0.1;

const double yb = -3.0;
const double ye = 3.0;
const double dy = 0.1;

const double zb = -3.0;
const double ze = 3.0;
const double dz = 0.1;

const double Tstart = 0.0;
const double Tmax = 40.0;
const double Tstep = 0.4;


const double err_abs = 1.0e-5;
const double err_rel = 1.0e-3;


//chua
const double C1_inv = 8;
const double L_inv = 7.0;
const double G = 0.7;
const double m0 = -0.5;
const double m1 = -0.8;

const double koef_mean = 0.7;

int num_set = 2;


typedef boost::array<double, 3> state_type;
vector<state_type> mas;
typedef runge_kutta_cash_karp54<state_type> error_stepper_rkck54;

vector<state_type> XYZgrid;

struct push_back_state_and_time {
    std::vector<state_type> &m_states;
    std::vector<double> &m_times;

    push_back_state_and_time(std::vector<state_type> &states, std::vector<double> &times)
            : m_states(states), m_times(times) {}

    void operator()(const state_type &x, double t) {
        m_states.push_back(x);
        m_times.push_back(t);
    }
};

void lorenz(const state_type &x, state_type &dxdt, double t) {
    dxdt[0] = sigma * (x[1] - x[0]);
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = -b * x[2] + x[0] * x[1];
}

double g(double x, double m0, double m1) {
    if ((x >= -1) and (x <= 1)) {
        return m1 * x;
    }
    if (x > 1) {
        return m0 * x + (m1 - m0);
    }
    if (x < -1) {
        return m0 * x - (m1 - m0);
    }
}

void chua(const state_type &x, state_type &dxdt, double t) {
    dxdt[0] = C1_inv * G * (x[1] - x[0]) - g(x[0], m0, m1);
    dxdt[1] = G * (x[0] - x[1]) + x[2];
    dxdt[2] = -L_inv * x[1];
}


void write_lorenz(const state_type &x, const double t) {
    ofstream fout;
    fout.open("test.txt", ios_base::app);
    cout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
    fout << t << '\t' << x[0] << '\t' << x[1] << '\t' << x[2] << endl;
    fout.close();
}

void write_lorenz2(const state_type &x, const double t) {
    mas.push_back({x[0], x[1], x[2]});
}

void write_data(vector<vector<state_type> > data) {
    FILE *fout = fopen("data1.dat", "wb");
    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[i].size(); j++) {
            fwrite(&data[i][j][0], sizeof(double), 1, fout);
            fwrite(&data[i][j][1], sizeof(double), 1, fout);
            fwrite(&data[i][j][2], sizeof(double), 1, fout);
        }
    }
    fclose(fout);
}

void write_data_mat(vector<vector<state_type> > data, int num_set) {

    char name[20];
    sprintf(name, "CHUA_%d.mat", num_set);

    std::cout << data[0].size() << " " << data[0][0].size() << " " << data.size() << endl;

    ofstream fout(name);

    fout << "# name: TrajBundle" << endl;
    fout << "# type: matrix" << endl;
    fout << "# ndims: 3" << endl;
    fout << data[0].size() << " " << data[0][0].size() << " " << data.size() << endl;

    for (int i = 0; i < data.size(); i++) {
        for (int j = 0; j < data[0][0].size(); j++) {
            for (int t = 0; t < data[0].size(); t++) {
                fout << data[i][t][j] << endl;
            }
        }
    }

    fout.close();
}


void getTrayDensity(vector<vector<state_type> > data) {

    double xbb, xee, ybb, yee, zbb, zee;
    vector<unsigned long int> outdata(XYZgrid.size());

    for (unsigned long long int i = 0; i < XYZgrid.size(); i++) {
        xbb = XYZgrid[i][0] - dx / 2;
        xee = XYZgrid[i][0] + dx / 2;

        ybb = XYZgrid[i][1] - dy / 2;
        yee = XYZgrid[i][1] + dy / 2;

        zbb = XYZgrid[i][2] - dz / 2;
        zee = XYZgrid[i][2] + dz / 2;
        if (i == 0) {
            cout << xbb << " " << XYZgrid[i][0] << " " << xee << endl;
            cout << ybb << " " << XYZgrid[i][1] << " " << yee << endl;
            cout << zbb << " " << XYZgrid[i][2] << " " << zee << endl;

        }

        for (unsigned long long int j = 0; j < data.size(); j++) {
            for (unsigned long int t = 0; t < data[0].size(); t++) {

                double xx = data[j][t][0];
                double yy = data[j][t][1];
                double zz = data[j][t][2];

                if ((xx >= xbb and xx <= xee) and
                    (yy >= ybb and yy <= yee) and
                    (zz >= zbb and zz <= zee)) {
                    outdata[i] += 1;
                }
            }

        }

    }

    for (unsigned long int i = 0; i < 10; i++) {
        cout << outdata[i] << endl;
    }

}

void getTrayDensityI(vector<vector<state_type> > data, int nx, int ny, int nz, vector<unsigned long int> &outdata2) {
    long long int ix, iy, iz;
    outdata2.resize(data.size());
    int count = 0;
    for (unsigned long long int j = 0; j < data.size(); j++) {
        for (unsigned long int t = 0; t < data[0].size(); t++) {

            ix = floor((data[j][t][0] - xb) / dx);
            iy = floor((data[j][t][1] - yb) / dy);
            iz = floor((data[j][t][2] - zb) / dz);
            if ((ix >= 0 and ix < nx) and (iy >= 0 and iy < ny) and (iz >= 0 and iz < nz)) {
                outdata2[ix + iy * nx + iz * nx * ny] += 1;
            }
        }
    }
}


void densFilter(vector<unsigned long int> &data, unsigned long int M) {

    for (int i = 0; i < data.size(); i++) {
        if (data[i] < M) {
            data[i] = 0;
        }
    }
}

unsigned long int meanM(vector<unsigned long int> data) {
    unsigned long int sum = 0;
    for (int i = 0; i < data.size(); i++) {
        sum += data[i];
    }
    return (unsigned long int) sum / data.size();
}

void write_data_dens_mat(vector<unsigned long int> data, int num_set) {

    char name[20];
    sprintf(name, "ITER_1_dens_chua_CHUA_%d.mat", num_set);

    ofstream fout(name);
    fout << "# name: ITER1" << endl;
    fout << "# type: matrix" << endl;
    fout << "# rows:" << data.size() << endl;
    fout << "# colums:" << 1 << endl;

    //fout<<data.size()<<endl;
    for (int i = 0; i < data.size(); i++) {
        fout << data[i] << endl;
    }
    fout.close();
}

void write_data_gridRef_mat(vector<state_type> data, int num_set) {
//	char *name;

    char name[20];
    sprintf(name, "CHUAref_%d.mat", num_set);

    ofstream fout(name);
    fout << "# name: GridRef" << endl;
    fout << "# type: matrix" << endl;
    fout << "# rows: " << data.size() << endl;
    fout << "# colums: " << data[0].size() << endl;
    for (int i = 0; i < data.size(); i++) {
        fout << data[i][0] << " " << data[i][1] << " " << data[i][2] << endl;
    }
    fout.close();
}


void filterOfGridXYZ(vector<unsigned long int> dens, int Nx, int Ny, int Nz, int num_set) {
    vector<state_type> XYZgrid_ref;
    unsigned long int num_cell = 0;
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                if (dens[num_cell] != 0) {
                    XYZgrid_ref.push_back(
                            {(xb + dx / 2) + x * (dx), (yb + dy / 2) + y * (dy), (zb + dz / 2) + z * (dz)});
                }
                num_cell++;
            }
        }
    }

    write_data_gridRef_mat(XYZgrid_ref, num_set);
}

//-------------------------------------------------------------------------------------------------------------------------
int main(int argc, char **argv) {
    vector<state_type> XYZgrid_s;
    int Nx = (xe - xb) / dx;
    int Ny = (ye - yb) / dy;
    int Nz = (ze - zb) / dz;
    for (int z = 0; z < Nz; z++) {
        for (int y = 0; y < Ny; y++) {
            for (int x = 0; x < Nx; x++) {
                XYZgrid_s.push_back({(xb + dx / 2) + x * (dx), (yb + dy / 2) + y * (dy), (zb + dz / 2) + z * (dz)});
            }
        }
    }
    XYZgrid = XYZgrid_s;

    cout << Nx << " " << Ny << " " << Nz << " " << XYZgrid.size() << endl;

    for (int i = 0; i < 10; i++) {
        cout << i << " " << XYZgrid[i][0] << "\t" << XYZgrid[i][1] << "\t" << XYZgrid[i][2] << endl;
    }
    size_t steps;
    vector<vector<state_type> > out_data(XYZgrid.size());
    int nt = omp_get_max_threads();
    //omp_set_num_threads(nt);
    cout << "omp = " << nt << endl;

    unsigned long long int i;
    size_t tstart = time(NULL);
    vector<state_type> x_vec;
    vector<double> times;
#pragma omp parallel shared(out_data, XYZgrid) private(i, x_vec, times)
    {
#pragma omp for schedule(static)

        for (i = 0; i < XYZgrid.size(); i++) {
            x_vec.clear();
            //	steps = integrate_const(make_controlled( err_abs , err_rel , error_stepper_rkck54() ), lorenz , XYZgrid[i] , Tstart , Tmax , Tstep , push_back_state_and_time( x_vec , times ));
            steps = integrate_const(make_controlled(err_abs, err_rel, error_stepper_rkck54()), chua, XYZgrid_s[i],
                                    Tstart, Tmax, Tstep, push_back_state_and_time(x_vec, times));
            out_data[i] = x_vec;
        }
    }
    cout << "Time = " << time(NULL) - tstart << endl;
    cout << "size_out_data = " << out_data.size() << endl;
    std::cout << out_data[0].size() << " " << out_data[0][0].size() << " " << out_data.size() << endl;

    write_data_mat(out_data, num_set);

    cout << "FINISH!!!" << endl;


//getTrayDensity(out_data);
    vector<unsigned long int> outdata;
    getTrayDensityI(out_data, Nx, Ny, Nz, outdata);
    for (unsigned long int i = 0; i < 10; i++) {
        cout << outdata[i] << endl;
    }
    unsigned long int M = meanM(outdata);
    cout << M << endl;

    densFilter(outdata, koef_mean * M);
    write_data_dens_mat(outdata, num_set);
    filterOfGridXYZ(outdata, Nx, Ny, Nz, num_set);

//	system("PAUSE");
}
