#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <vector>
#include <omp.h>
#include <ctime>
#include <string>
#include <cmath>

#include <unordered_map>
#include <iostream>

namespace std
{
    template <>
    struct hash<std::pair<size_t, size_t>>
    {
        size_t operator() (const std::pair<size_t, size_t>& val) const
        {
            return (val.first << ((sizeof(size_t) / 2) * 8)) + val.second;
        }
    };
}

using namespace std;
using namespace boost::numeric::odeint;

// space parameters:
//left and right ends of vector_matrix interval; vector_matrix step
const double xb = -5.0;
const double xe = 5.0;
const double dx = 0.1;

//left and right ends of y interval; y step
const double yb = -3.0;
const double ye = 3.0;
const double dy = 0.1;

//left and right ends of z interval; z step
const double zb = -3.0;
const double ze = 3.0;
const double dz = 0.1;

// time parameters: beginninig and end of time interval; time step
const double TStart = 0.0;
const double Tmax = 40.0;
const double Tstep = 1;


// maximum absolut and relative errors for nimerical solver of ODE
const double err_abs = 1.0e-5;
const double err_rel = 1.0e-3;


//right hand side of ODE system: chua circuit
const double C1_inv = 8;
const double L_inv = 7.0;
const double G = 0.7;
const double m0 = -0.5;
const double m1 = -0.8;

double g(double x, double m0, double m1) {
    if ((x >= -1) and (x <= 1)) {
        return m1 * x;
    }
    else if (x > 1) {
        return m0 * x + (m1 - m0);
    }
    else /* (x < -1) */ {
        return m0 * x - (m1 - m0);
    }
}

int num_set = 2;

//define type for an element of the state space - 3D vector
using state_type = std::array<double, 3>;
vector<state_type> mas;

using vector_matrix = std::unordered_map<std::pair<size_t, size_t>, state_type>;

void chua(const state_type &x, state_type &dxdt, double t) {
    dxdt[0] = C1_inv * G * (x[1] - x[0]) - g(x[0], m0, m1);
    dxdt[1] = G * (x[0] - x[1]) + x[2];
    dxdt[2] = -L_inv * x[1];
}




// ODE solver (see help for odeint, I dont know, what does this mean)
typedef runge_kutta_cash_karp54<state_type> error_stepper_rkck54;

// Space grid: coordinates of the centers of the rectangles.
// the center of the rectangle with the indeces i_x, j_y, l_z has the position
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

void write_data_mat(vector_matrix data, int num_set, size_t rows, size_t columns) {

    char name[20];
    sprintf(name, "CHUA_%d.mat", num_set);

//    std::cout << data[0].size() << " " << data[0][0].size() << " " << data.size() << endl;

    ofstream fout(name);

    fout << "# name: TrajBundle" << endl;
    fout << "# type: matrix" << endl;
    fout << "# ndims: 3" << endl;
//    fout << data[0].size() << " " << data[0][0].size() << " " << data.size() << endl;

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < 3; j++) {
            for (int t = 0; t < columns; t++) {
                fout << data[{i, t}][j] << endl;
            }
        }
    }

    fout.close();
}


void getTrayDensityT2T(vector_matrix data, int nx, int ny, int nz, vector<double> &DensT,
                       vector<double> &Dens2T, size_t rows, size_t columns) {
// input data - trajectory bundle
// input nx, ny, nz - sizes of initial data grid
// output DensT - density of trajectories on [0, T]
// output DensT2 - density of trajectories on [0, 2T]
// DensT, DensT2 are one-dimesional arrays; density of trajectories in the rectangle with the index (ix,iy,iz)  has the position Dens2T[ix+iy*nx+iz*nx*ny]

    long long int ix, iy, iz;
    DensT.resize(rows);
    Dens2T.resize(rows);
    int size2T, sizeT, TrajNum;
    size2T = columns;
    sizeT = size2T / 2;
    TrajNum = rows;

    int count = 0;
    for (unsigned long long int j = 0; j < TrajNum; j++) {
        for (unsigned long int t = 0; t < sizeT; t++) {
            ix = floor((data[{j, t}][0] - xb) / dx);
            iy = floor((data[{j, t}][1] - yb) / dy);
            iz = floor((data[{j, t}][2] - zb) / dz);
            if ((ix >= 0 and ix < nx) and (iy >= 0 and iy < ny) and (iz >= 0 and iz < nz)) {
                DensT[ix + iy * nx + iz * nx * ny] += 1;
            }
        }
        for (unsigned long int t = sizeT; t < sizeT; t++) {
            ix = floor((data[{j, t}][0] - xb) / dx);
            iy = floor((data[{j, t}][1] - yb) / dy);
            iz = floor((data[{j, t}][2] - zb) / dz);
            if ((ix >= 0 and ix < nx) and (iy >= 0 and iy < ny) and (iz >= 0 and iz < nz)) {
                Dens2T[ix + iy * nx + iz * nx * ny] += 1;
            }
        }
    }
    for (unsigned long int i = 0; i < Dens2T.size(); i++) {
        Dens2T[i] += DensT[i];
    }
    for (unsigned long int i = 0; i < Dens2T.size(); i++) {
        DensT[i] = DensT[i] / (sizeT * TrajNum);
    }
    for (unsigned long int i = 0; i < Dens2T.size(); i++) {
        Dens2T[i] = Dens2T[i] / (size2T * TrajNum);
    }
}


void densFilterT2T(vector<double> &dataT, vector<double> &data2T) {
// input dataT - density of trajectories on [0, T]
// input data2T - density of trajectories on [0, 2T]
// The result of filtration is written to data2T
    for (unsigned long int i = 0; i < dataT.size(); i++) {
        if (data2T[i] < 0.8 * dataT[i]) {
            data2T[i] = 0;
        }
    }
}

[i, j, k]
{1, 0, 2}
1 1 0 = 1
1 1 2 = 2
void write_data_dens_mat(vector<double> data, int num_set) {
//write down a matrix to a file in the octave readable format
    char name[35];
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
// writes a space grid to a file
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


void filterOfGridXYZ(vector<double> dens, int Nx, int Ny, int Nz, int num_set) {
//not tested!
// filters grid accordin to density of trajectories matrix and write it to a file
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
int main(int argc, char **argv)
{

    //set of initial data for trajectory bundle
    vector<state_type> XYZgrid_s;

    //calculation of the initial data grid
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

    //trajectory bundle
//    vector<vector<state_type> > out_data(XYZgrid.size());
    vector_matrix out_data;

    int nt = omp_get_max_threads();
    //omp_set_num_threads(nt);
    cout << "omp = " << nt << endl;

    size_t rows_count = XYZgrid.size();
    size_t columns_count;

    unsigned long long int i;
    size_t tstart = time(NULL);
    vector<state_type> x_vec;
    vector<double> times;

#pragma omp parallel shared(out_data, XYZgrid) private(i, x_vec, times)
    {
#pragma omp for schedule(static)

        for (i = 0; i < XYZgrid.size(); i++) {
            x_vec.clear();
            //	steps = integrate_const(make_controlled( err_abs , err_rel , error_stepper_rkck54() ), lorenz , XYZgrid[i] , TStart , Tmax , Tstep , push_back_state_and_time( x_vec , times ));
            steps = integrate_const(make_controlled(err_abs, err_rel, error_stepper_rkck54()), chua, XYZgrid_s[i],
                                    TStart, Tmax, Tstep, push_back_state_and_time(x_vec, times));

            columns_count = x_vec.size();

            for (int j = 0; j < columns_count; ++j)
            {
                out_data.insert({{i, j}, x_vec[j]});
//                out_data[{i, j}] = x_vec[j];
            }
        }
    }
    cout << "Time = " << time(NULL) - tstart << endl;
    cout << "size_out_data = " << out_data.size() << endl;
//    std::cout << out_data[0].size() << " " << out_data[0][0].size() << " " << out_data.size() << endl;

    write_data_mat(out_data, num_set, rows_count, columns_count);

    cout << "FINISH!!!" << endl;


    // this part is not tested!
    //getTrayDensity(out_data);
    vector<double> dens_data_T;
    vector<double> dens_data_2T;
    getTrayDensityT2T(out_data, Nx, Ny, Nz, dens_data_T, dens_data_2T, rows_count, columns_count);


    densFilterT2T(dens_data_T, dens_data_2T);
    write_data_dens_mat(dens_data_2T, num_set);
    filterOfGridXYZ(dens_data_2T, Nx, Ny, Nz, num_set);


    return 0;
//	system("PAUSE");
}
