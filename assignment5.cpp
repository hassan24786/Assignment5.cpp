/*Description: This code creates multiple particle objects, derived from a single lepton class. 
Each particle has its own momentum and fundamental properties, which can be printed. 
The code also involves basic interaction of particles with a calorimeter. 
Author: Hassan Hashmi
Date: 20/04/2024*/ 

#include <iostream>
#include <string>
#include <iomanip>
#include <vector>
#include <memory> 

const double speed_of_light = 1; // in natural units

class FourMomentum 
{
private:                                                    
  std::unique_ptr<std::vector<double>>momentum_vector;     

public:                                                      
  FourMomentum(double input_energy, double input_momentum_x, double input_momentum_y, double input_momentum_z):
  momentum_vector{std::make_unique<std::vector<double>>(4)}
  {
    set_energy(input_energy);                               
    set_momentum_x(input_momentum_x);
    set_momentum_y(input_momentum_y);
    set_momentum_z(input_momentum_z);
  }
  void set_energy(double input_energy)                                  
  {
  if(input_energy < 0)
    std::cout<<"The energy value cannot be negative."<<std::endl;
  else
    (*momentum_vector)[0] = input_energy / speed_of_light;
  }
  void set_momentum_x(double input_momentum_x) {(*momentum_vector)[1] = input_momentum_x;}
  void set_momentum_y(double input_momentum_y) {(*momentum_vector)[2] = input_momentum_y;}
  void set_momentum_z(double input_momentum_z) {(*momentum_vector)[3] = input_momentum_z;}
  double get_energy() const {return(*momentum_vector)[0];}            
  double get_momentum_x() const {return (*momentum_vector)[1];}
  double get_momentum_y() const {return (*momentum_vector)[2];}
  double get_momentum_z() const {return (*momentum_vector)[3];}
  FourMomentum(const FourMomentum &particle_to_copy)                 
  {
    momentum_vector = std::make_unique<std::vector<double>>(*particle_to_copy.momentum_vector);
  }                             
  FourMomentum &operator=(const FourMomentum &particle_to_copy)     
  {
    std::cout<<"\nCalling momentum vector copy assignment operator. "<<std::endl;
    if(&particle_to_copy == this) 
      return *this;
    momentum_vector = std::make_unique<std::vector<double>>(*particle_to_copy.momentum_vector);
    return *this;
  }
  FourMomentum(FourMomentum &&particle_to_move)       
  {
    std::cout<<"\nCalling momentum vector move constructor. "<<std::endl;
    momentum_vector = std::move(particle_to_move.momentum_vector);
  }            
  FourMomentum &operator=(FourMomentum &&particle_to_move) 
  {
    std::cout<<"\nCalling momentum vector move assignment operator. "<<std::endl;
    momentum_vector = std::move(particle_to_move.momentum_vector);
    return *this; 
  }
  ~FourMomentum() {}
};

class lepton 
{
protected:         
  double rest_mass;
  int charge;
  std::unique_ptr<FourMomentum> momentum; 
public:
  lepton() = default; 
  lepton(double input_mass, int input_charge, double input_energy, double input_momentum_x, 
  double input_momentum_y, double input_momentum_z):
  rest_mass{input_mass}, charge{input_charge}, momentum{std::make_unique<FourMomentum>(input_energy, 
  input_momentum_x, input_momentum_y, input_momentum_z)} {} 
  const FourMomentum &get_momentum() const {return *momentum;} 
  friend std::vector<double> sum_of_momentum(const FourMomentum &particle_one_momentum, 
  const FourMomentum &particle_two_momentum);
  friend double momentum_dot_product(const FourMomentum &particle_one_momentum, const FourMomentum &particle_two_momentum);
  virtual ~lepton() = default;
  lepton(const lepton &particle_to_copy) 
  {
    momentum = std::make_unique<FourMomentum>(*particle_to_copy.momentum);
  }
  lepton &operator=(const lepton &particle_to_copy) 
  {
    std::cout<<"\nCalling lepton copy assignment operator. "<<std::endl;
    if(&particle_to_copy== this) 
      return *this;
    momentum = std::make_unique<FourMomentum>(*particle_to_copy.momentum);
    return *this;
  }
  lepton(lepton &&particle_to_move)
  {
    std::cout<<"\nCalling lepton move constructor. "<<std::endl;
    rest_mass=particle_to_move.rest_mass;
    charge = particle_to_move.charge ;
    momentum = std::move(particle_to_move.momentum);
    particle_to_move.rest_mass = 0;
    particle_to_move.charge = 0;
  }
  lepton &operator=(lepton &&particle_to_move) 
  {
    std::cout<<"\nCalling lepton move assignment operator. Particle properties after move: "<<std::endl;
    rest_mass = particle_to_move.rest_mass;
    charge = particle_to_move.charge;
    momentum = std::move(particle_to_move.momentum); 
    particle_to_move.rest_mass = 0;
    particle_to_move.charge = 0;
    return *this;
  }
  friend void print_lepton_data(const lepton &input_lepton);
};

class Calorimeter 
{
private:
  double EM1, EM2, HAD1, HAD2;  //calorimeter layers
public:
  Calorimeter() = default;
  void set_calorimeter_energies(const std::vector<double> &calorimeter_energies) 
  {
    EM1 = calorimeter_energies[0];
    EM2 = calorimeter_energies[1];
    HAD1 = calorimeter_energies[2];
    HAD2 = calorimeter_energies[3];
  }
  double total_calorimeter_energy() const {return EM1 + EM2 + HAD1 + HAD2;}
};

class Electron : public lepton 
{
private:
  Calorimeter calorimeter;
  std::vector<double> lepton_energies;
public:
  Electron(double input_energy, double input_momentum_x, double input_momentum_y, double input_momentum_z,
  const std::vector<double>& energies):
  lepton(0.511, -1, input_energy, input_momentum_x, input_momentum_y, input_momentum_z), lepton_energies{energies}
  {
    calorimeter.set_calorimeter_energies(lepton_energies);
    if(std::abs(get_momentum().get_energy() - calorimeter.total_calorimeter_energy()) > 1e-4)
      std::cout<<"The total energy deposited by the electron in the calorimeter does not match the electron's four momentum energy. "<<std::endl; 
  }
  std::vector<double> get_lepton_energies() const {return lepton_energies;}
  friend void print_electron_data(const Electron &electron_particle);
};

class Muon : public lepton 
{
private:
  bool muon_isolation;
public:
  Muon(double input_energy, double input_momentum_x, double input_momentum_y, double input_momentum_z, bool isolation_status):
  lepton(105.7, -1, input_energy, input_momentum_x, input_momentum_y, input_momentum_z), muon_isolation{isolation_status} {}
  bool get_isolation_status() const {return muon_isolation;}
  void set_isolation_status (bool isolation_status) {muon_isolation = isolation_status;}
  friend void print_muon_data(const Muon &muon_particle);
};

class Neutrino : public lepton 
{
private:
  bool hasInteracted;
  std::string flavour; 
public:
  Neutrino(double input_energy, double input_momentum_x, double input_momentum_y, double input_momentum_z, 
  std::string input_flavour, bool interaction_status):
  lepton(0, 0, input_energy, input_momentum_x, input_momentum_y, input_momentum_z), flavour{input_flavour}, hasInteracted{interaction_status} {}
  void set_interaction_status(bool interaction_status) {hasInteracted = interaction_status;}
  bool get_interaction_status() const {return hasInteracted;}
  const std::string &get_flavour() const {return flavour;}
  friend void print_neutrino_data(const Neutrino &neutron_particle);
};

class Tau : public lepton 
{
private:
  bool leptonic_decay; 
  std::vector<std::shared_ptr<lepton>>decay_products; 
public:
  Tau(double input_energy, double input_momentum_x, double input_momentum_y, double input_momentum_z, bool leptonic_or_hadronic):
  lepton(1777, -1, input_energy, input_momentum_x, input_momentum_y, input_momentum_z), leptonic_decay{leptonic_or_hadronic} {}
  void add_decay_product(const std::shared_ptr<lepton> &product) {decay_products.push_back(product);}
  bool get_is_leptonic() const {return leptonic_decay;}
  friend void print_tau_data(const Tau &tau_particle);
};

void print_tau_data(const Tau &tau_particle) 
{
  std::cout<<"Tau particle information:"<<std::endl;
  print_lepton_data(tau_particle);
  std::cout<<"\nCurrent decay mode: "<<std::endl;
  if(tau_particle.get_is_leptonic ()) 
    std::cout<<"The tau particle decays leptonically. "<<std::endl;
  else 
    std::cout<<"The tau particle decays hadronically. "<<std::endl;
  std::cout<<"The decay produces the following particles :"<<std::endl;
  for(const std::shared_ptr<lepton> &decay_particle : tau_particle.decay_products) 
  {
    print_lepton_data(*decay_particle);
  }
}

void print_neutrino_data(const Neutrino &neutrino_particle) 
{
  std::cout<<"Neutrino particle information: "<<std::endl;
  print_lepton_data(neutrino_particle); 
  std::cout<<"Flavour of the neutrino: "<<neutrino_particle.get_flavour()<<std::endl;
  std::cout<<"Interaction status (true = isolated, false = not isolated): "<<neutrino_particle.get_interaction_status()<<std::endl;
}

void print_electron_data(const Electron &electron_particle) 
{
  print_lepton_data(electron_particle); 
  std::cout<<"Deposited energies in each layer of the calorimeter: ";
  for(double energy : electron_particle.get_lepton_energies()) 
    {std::cout<<energy<<"MeV ";}
}

void print_muon_data(const Muon &muon_particle) 
{
  std::cout<<"Muon particle information: "<<std::endl;
  print_lepton_data(muon_particle); 
  std::cout<<"Isolation status (true=isolated, false = not isolated): "<<muon_particle.get_isolation_status()<<std::endl;
}

void print_lepton_data(const lepton &input_lepton)      
{
  std::cout<<std::fixed<<std::setprecision(3);
  std::cout<<"Rest Mass (MeV): "<<input_lepton.rest_mass<<std::endl;
  std::cout<<"Charge: "<<input_lepton.charge<<std::endl;
  const FourMomentum &momentum = input_lepton.get_momentum();
  std::cout<<"E/c: "<<momentum.get_energy()<<std::endl;
  std::cout<<"Momentum in x direction: "<<momentum.get_momentum_x()<<std::endl;
  std::cout<<"Momentum in y direction: "<<momentum.get_momentum_y()<<std::endl;
  std::cout<<"Momentum in z direction: "<<momentum.get_momentum_z()<<std::endl;
}

std::vector<double> sum_of_momentum(const FourMomentum &particle_one_momentum, const FourMomentum &particle_two_momentum) 
{
  std::vector<double> sum_vector(4, 0.0);
  sum_vector[0] = particle_one_momentum.get_energy() + particle_two_momentum.get_energy();
  sum_vector[1] = particle_one_momentum.get_momentum_x() + particle_two_momentum.get_momentum_x();
  sum_vector[2] = particle_one_momentum.get_momentum_y() + particle_two_momentum.get_momentum_y();
  sum_vector[3] = particle_one_momentum.get_momentum_z() + particle_two_momentum.get_momentum_z();
  return sum_vector;
}

double momentum_dot_product(const FourMomentum &particle_one_momentum, const FourMomentum &particle_two_momentum) 
{
  double dot_product = particle_one_momentum.get_energy() * particle_two_momentum.get_energy() +
                       particle_one_momentum.get_momentum_x() * particle_two_momentum.get_momentum_x() +
                       particle_one_momentum.get_momentum_y() * particle_two_momentum.get_momentum_y() +
                       particle_one_momentum.get_momentum_z() * particle_two_momentum.get_momentum_z();
  return dot_product;
}

int main() 
{
  std::vector<lepton> lepton_information_vector = 
  {
    lepton(0.511, -1, -1, 0, 1, 0),
    lepton(0.511, -1, 0, 2, 0, 0),
    lepton(105.7, -1, 0, 0, 0, 0),
    lepton(105.7, -1, 0, 0, 0, 0),
    lepton(105.7, -1, 0, 0, 0, 0),
    lepton(105.7, -1, 0, 0, 0, 0),
    lepton(0.511, 1, 0, 0, 1, 0),
    lepton(105.7, 1, 0, 0, 1, 0),
    lepton(0, 0, 0, 0, 0, 0), 
    lepton(0, 0, 0, 0, 0, 0) 
  };
  std::vector<double> electron_sum = sum_of_momentum(lepton_information_vector[0].get_momentum(), lepton_information_vector[1].get_momentum());
  std::cout<<"Sum of electron four momentum vectors - "<<std::endl;
  std::cout<<"E/c : "<<electron_sum[0]<<std::endl;
  std::cout<<"Momentum in x direction: "<<electron_sum[1]<<std::endl; 
  std::cout<<"Momentum in y direction: "<<electron_sum[2]<<std::endl;  
  std::cout<<"Momentum in z direction: "<<electron_sum[3]<<std::endl;   
  std::cout<<std::endl;
  double dot_product = momentum_dot_product(lepton_information_vector[6].get_momentum(), lepton_information_vector[7].get_momentum());
  std::cout<<"Dot product of the antielectron and antimuon four momentum vectors: "<<dot_product << std::endl;
  std::vector<double> lepton_energies = {1,2,3,4}; //example energy values (just for test)
  std::unique_ptr<Electron> new_electron_pointer{new Electron(50, 3, 2, 0, lepton_energies)};
  Electron electron_two(0.511, 0, 0, 0, lepton_energies); 
  electron_two = std::move(*new_electron_pointer);
  print_electron_data(electron_two);
  std::shared_ptr<Tau> tau_pointer {new Tau(2, 0.5, 0, 0, true)}; 
  std::shared_ptr<lepton> detector1_tau;
  std::shared_ptr<lepton> detector2_tau;
  detector1_tau = tau_pointer;
  detector2_tau = tau_pointer;
  return 0; 
}