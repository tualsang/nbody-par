#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

#include <omp.h>    // Open

double G = 6.674*std::pow(10,-11);
//double G = 1;

struct simulation {
  std::size_t nbpart;    // Number of particles.
  
  std::vector<double> mass;   // Mass of Particles

  //position
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;

  //velocity
  std::vector<double> vx;
  std::vector<double> vy;
  std::vector<double> vz;

  //force
  std::vector<double> fx;
  std::vector<double> fy;
  std::vector<double> fz;

  
  simulation(std::size_t nb)
    :nbpart(nb), mass(nb),
     x(nb,0.0), y(nb), z(nb),
     vx(nb), vy(nb), vz(nb),
     fx(nb), fy(nb), fz(nb) 
  {}
};


void random_init(simulation& s) {
  std::random_device rd;  
  std::mt19937 gen(rd());
  std::uniform_real_distribution dismass(0.9, 1.);
  std::normal_distribution dispos(0.0, 1.0);
  std::normal_distribution disvel(0.0, 1.0);

  for (std::size_t i = 0; i<s.nbpart; ++i) {
    s.mass[i] = dismass(gen);

    s.x[i] = dispos(gen);
    s.y[i] = dispos(gen);
    s.z[i] = dispos(gen);
    s.z[i] = 0.;
    
    s.vx[i] = disvel(gen);
    s.vy[i] = disvel(gen);
    s.vz[i] = disvel(gen);
    s.vz[i] = 0.;
    s.vx[i] = s.y[i]*1.5;
    s.vy[i] = -s.x[i]*1.5;
  }

  return;
  //normalize velocity (using normalization found on some physicis blog)
  double meanmass = 0;
  double meanmassvx = 0;
  double meanmassvy = 0;
  double meanmassvz = 0;
  for (std::size_t i = 0; i<s.nbpart; ++i) {
    meanmass += s.mass[i];
    meanmassvx += s.mass[i] * s.vx[i];
    meanmassvy += s.mass[i] * s.vy[i];
    meanmassvz += s.mass[i] * s.vz[i];
  }
  for (std::size_t i = 0; i<s.nbpart; ++i) {
    s.vx[i] -= meanmassvx/meanmass;
    s.vy[i] -= meanmassvy/meanmass;
    s.vz[i] -= meanmassvz/meanmass;
  }
  
}

void init_solar(simulation& s) {
  enum Planets {SUN, MERCURY, VENUS, EARTH, MARS, JUPITER, SATURN, URANUS, NEPTUNE, MOON};
  s = simulation(10);

  // Masses in kg
  s.mass[SUN] = 1.9891 * std::pow(10, 30);
  s.mass[MERCURY] = 3.285 * std::pow(10, 23);
  s.mass[VENUS] = 4.867 * std::pow(10, 24);
  s.mass[EARTH] = 5.972 * std::pow(10, 24);
  s.mass[MARS] = 6.39 * std::pow(10, 23);
  s.mass[JUPITER] = 1.898 * std::pow(10, 27);
  s.mass[SATURN] = 5.683 * std::pow(10, 26);
  s.mass[URANUS] = 8.681 * std::pow(10, 25);
  s.mass[NEPTUNE] = 1.024 * std::pow(10, 26);
  s.mass[MOON] = 7.342 * std::pow(10, 22);

  // Positions (in meters) and velocities (in m/s)
  double AU = 1.496 * std::pow(10, 11); // Astronomical Unit

  s.x = {0, 0.39*AU, 0.72*AU, 1.0*AU, 1.52*AU, 5.20*AU, 9.58*AU, 19.22*AU, 30.05*AU, 1.0*AU + 3.844*std::pow(10, 8)};
  s.y = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  s.z = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

  s.vx = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  s.vy = {0, 47870, 35020, 29780, 24130, 13070, 9680, 6800, 5430, 29780 + 1022};
  s.vz = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
}

//meant to update the force that from applies on to
void update_force(simulation& s, std::size_t from, std::size_t to) {
  double softening = .1;
  double dist_sq = std::pow(s.x[from]-s.x[to],2)
    + std::pow(s.y[from]-s.y[to],2)
    + std::pow(s.z[from]-s.z[to],2);
  double F = G * s.mass[from]*s.mass[to]/(dist_sq+softening); //that the strength of the force

  //direction
  double dx = s.x[from]-s.x[to];
  double dy = s.y[from]-s.y[to];
  double dz = s.z[from]-s.z[to];
  double norm = std::sqrt(dx*dx+dy*dy+dz*dz);
  
  dx = dx/norm;
  dy = dy/norm;
  dz = dz/norm;

  //apply force
  s.fx[to] += dx*F;
  s.fy[to] += dy*F;
  s.fz[to] += dz*F;
}

void compute_forces(simulation& s) {
  #pragma omp parallel
  {
    std::vector<double> fx(s.nbpart, 0.), fy(s.nbpart, 0.), fz(s.nbpart, 0.);

    #pragma omp for
    for (std::size_t i=0; i<s.nbpart; ++i) {
      for (std::size_t j = 0; j<s.nbpart; ++j) {
        if (i != j) {
          double softening = .1;
          double dist_sq = std::pow(s.x[i]-s.x[j],2)
          + std::pow(s.y[i]-s.y[j], 2)
          + std::pow(s.z[i]-s.z[j], 2);
          double F= G * s.mass[i]*s.mass[j]/(dist_sq+softening);

          // direction
          double dx = s.x[i] - s.x[j];
          double dy = s.y[i] - s.y[j];
          double dz = s.z[i] - s.z[j];
          double norm = std::sqrt(dx*dx+dy*dy+dz*dz);

          dx = dx/norm;
          dy = dy/norm;
          dz = dz/norm;

          //apply force
          fx[i] += dx*F;
          fy[i] += dy*F;
          fz[i] += dz*F;
        }
      }
    }
    #pragma omp critical
    {
      for (std::size_t i = 0; i<s.nbpart; ++i) {
        s.fx[i] += fx[i];
        s.fy[i] += fy[i];
        s.fz[i] += fz[i];
      }
    }
  }
}

void reset_force(simulation& s) {
  for (size_t i=0; i<s.nbpart; ++i) {
    s.fx[i] = 0.;
    s.fy[i] = 0.;
    s.fz[i] = 0.;
  }
}

void apply_force(simulation& s, std::size_t i, double dt) {
  s.vx[i] += s.fx[i]/s.mass[i]*dt;
  s.vy[i] += s.fy[i]/s.mass[i]*dt;
  s.vz[i] += s.fz[i]/s.mass[i]*dt;
}

void update_position(simulation& s, std::size_t i, double dt) {
  s.x[i] += s.vx[i]*dt;
  s.y[i] += s.vy[i]*dt;
  s.z[i] += s.vz[i]*dt;
}

void update_velocities(simulation& s, double dt) {
  #pragma omp parallel for
  for (std::size_t i = 0; i < s.nbpart; ++i) {
    apply_force(s, i, dt);
  }
}

void update_positions(simulation&s, double dt) {
  #pragma omp parallel for
  for (std::size_t i = 0; i < s.nbpart; ++i) {
    update_position(s, i, dt);
  }
}

// Honestly we could combine the update_velocities & update_positions into 1 method: update_particles.
void update_particles(simulation& s, double dt) {
  #pragma omp parallel for
  for (std::size_t i = 0; i < s.nbpart; ++i) {
    apply_force(s, i, dt);
    update_position(s, i, dt);
  }
}

void dump_state(simulation& s) {
  std::cout<<s.nbpart<<'\t';
  for (std::size_t i=0; i<s.nbpart; ++i) {
    std::cout<<s.mass[i]<<'\t';
    std::cout<<s.x[i]<<'\t'<<s.y[i]<<'\t'<<s.z[i]<<'\t';
    std::cout<<s.vx[i]<<'\t'<<s.vy[i]<<'\t'<<s.vz[i]<<'\t';
    std::cout<<s.fx[i]<<'\t'<<s.fy[i]<<'\t'<<s.fz[i]<<'\t';
  }
  std::cout<<'\n';
}

void load_from_file(simulation& s, std::string filename) {
  std::ifstream in (filename);
  std::size_t nbpart;
  in>>nbpart;
  s = simulation(nbpart);
  for (std::size_t i=0; i<s.nbpart; ++i) {
    in>>s.mass[i];
    in >>  s.x[i] >>  s.y[i] >>  s.z[i];
    in >> s.vx[i] >> s.vy[i] >> s.vz[i];
    in >> s.fx[i] >> s.fy[i] >> s.fz[i];
  }
  if (!in.good())
    throw "kaboom";
}

int main(int argc, char* argv[]) {
  if (argc != 5) {
    std::cerr
      <<"usage: "<<argv[0]<<" <input> <dt> <nbstep> <printevery>"<<"\n"
      <<"input can be:"<<"\n"
      <<"a number (random initialization)"<<"\n"
      <<"planet (initialize with solar system)"<<"\n"
      <<"a filename (load from file in singleline tsv)"<<"\n";
    return -1;
  }
  
  double dt = std::atof(argv[2]); //in seconds
  std::size_t nbstep = std::atol(argv[3]);
  std::size_t printevery = std::atol(argv[4]);
  
  simulation s(1);

  //parse command line
  {
    std::size_t nbpart = std::atol(argv[1]); //return 0 if not a number
    if ( nbpart > 0) {
      s = simulation(nbpart);
      random_init(s);
    } else {
        std::string inputparam = argv[1];
        if (inputparam == "planet") {
	      init_solar(s);
      } else  {
    	  load_from_file(s, inputparam);
      }
    }    
  }

  
  for (std::size_t step = 0; step< nbstep; step++) {
    if (step % printevery == 0)
      dump_state(s);
  
    reset_force(s);
  //   for (std::size_t i=0; i<s.nbpart; ++i)
  //     for (std::size_t j=0; j<s.nbpart; ++j)
	// if (i != j)
	//   update_force(s, i, j);
    compute_forces(s);

    // for (std::size_t i=0; i<s.nbpart; ++i) {
    //   apply_force(s, i, dt);
    //   update_position(s, i, dt);
    // }
    update_velocities(s, dt);
    update_positions(s, dt);
  }
  
  //dump_state(s);  
  return 0;
}
