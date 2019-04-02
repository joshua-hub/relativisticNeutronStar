// phys3071 as05 melsom 42593249

/******************************************************************************
Description: This program solves for the radius of a Neutron star using both the
Tolman-Oppenheimer-Volkoff (TOV) equation of state and the Newtonian equation 
of state. Te solutions are calculated over a range of central energy densities.

Inputs: The user is not required to input anything

Calculations: using a midpoint version of Euler's steps, the program integrates 
over the radius until the pressure of the star is less than zero. This is 
performed over both equations separately as pressure goes to zero at different 
rates. The internal energy is varied and the calculations are repeated.

Outputs: The program will create a file named "massvsrdius.dat" in the current
directory. The columns are initial energy density, calculated radius, 
calculated mass, calculated radius using TOV and calculated mass using TOV.

compiled as gcc as05-42593249-melsom.c -o ws09 -lm -Wall
******************************************************************************/

#include<math.h>
#include<stdio.h>
#include<stdlib.h>

// Define constants************************************************************
#define BAG 57
#define K 1.322e-42
#define SOLARMASS 1.98892e30
#define C 2.99e8
#define EV 1.602e-19

// Function prototypes*********************************************************
double dedr(double epsi, double mass, double r);
double dmdr(double epsi, double r);
double dedrTOV(double epsi, double mass, double r);

// Begin main function ********************************************************
int main () {
  int epsi_count;
  // dr time steps are this as my worksheet optimized by halving dr until 
  // radius didn't change much over halving. 3 iterations before it stopped.
  // all steps are divided into new, half new and old for the midpoint method
  // epsilon, mass, epsilon using TOV equation and mass using TOV epsilon
  // need to define a pressure and a solar mass for each method of calculation
  double r, rTOV, dr= (1.0e15/8.0), epsi_step;
  double epsi_old, epsi_new, epsi_half_new, epsi_zero;
  double epsiTOV_old, epsiTOV_new, epsiTOV_half_new;
  double mass_old, mass_new, mass_half_new;
  double massTOV_old, massTOV_new, massTOV_half_new;
  double pressure_new, mass_solar, pressureTOV_new, massTOV_solar;
  FILE *fp; // for pointing at file.
  fp = fopen("massvsradius.dat", "w");
  fprintf(fp,"#phys3071 as05 melsom 42593249\n");
  
  printf("\nThis program will create a data file named massvsradius.dat which\n"
        "contains calculated results for the radius of a neutron star for a\n"
        "range of central energy densities. Calculations will be done using \n"
        "Newtonian and a Relativistic formulas.\n");

  epsi_zero= 4.2;
  epsi_step= 0.2; 
  
                  // loops over epsilon = 4.2B to 20B  
                  // inside loop calculates solutions for both Newtonian and
                  // Relativistic solution of a Neutron star radius
  for (epsi_count= 0; epsi_count< 80; epsi_count++) {
  
    // clearing variables for Newtonian calculations.
    epsi_old= epsi_zero* BAG;
    mass_old= 0.0;
    pressure_new= 0.0;

               // solves radius for Newtonian neutron star
               // completes a half step and a full step.
    r= 1.0e16; // calculates pressure, radius, mass
    do {
      epsi_half_new= epsi_old+ (dr/2.0)* dedr(epsi_old, mass_old, r);
      mass_half_new= mass_old+ (dr/2.0)* dmdr(epsi_old, r);

      epsi_new= epsi_old+ dr* dedr(epsi_half_new, mass_half_new, r+dr);
      mass_new= mass_old+ dr* dmdr(epsi_half_new, r+dr);

      epsi_old= epsi_new;
      mass_old= mass_new;
      mass_solar = 1.0e6* EV* mass_new/ (pow(C,2.0)* SOLARMASS);
      pressure_new= (1.0/ 3.0)* (epsi_new- 4.0* BAG);
            
      r= r+ dr;
    } while (pressure_new> 0.0);

    // clearing variables for the TOV calculations
    epsiTOV_old= epsi_zero* BAG;
    massTOV_old= 0.0;
    pressureTOV_new= 0.0;

               // solves radius for relativity neutron star
               // completes a half step and a full step.
    rTOV= 1.0e16; // calculates pressure, radius, mass
    do {
      epsiTOV_half_new = epsiTOV_old+ (dr/2.0)* 
                         dedrTOV(epsiTOV_old, massTOV_old, rTOV);
      massTOV_half_new= massTOV_old+ (dr/2.0)* dmdr(epsiTOV_old, rTOV);
      
      epsiTOV_new = epsiTOV_old+ dr* 
                    dedrTOV(epsiTOV_half_new, massTOV_half_new, rTOV+dr);
      massTOV_new = massTOV_old+ dr* dmdr(epsiTOV_half_new, rTOV+dr);

      epsiTOV_old= epsiTOV_new;
      massTOV_old= massTOV_new;
      massTOV_solar = 1.0e6* EV* massTOV_new/ (pow(C,2.0)* SOLARMASS);   
      pressureTOV_new= (1.0/3.0)* (epsiTOV_new- 4.0* BAG);
      
      rTOV =rTOV +dr;
    } while (pressureTOV_new> 0.0);

  // program outputs to file "massvsradius.dat" the following columns.
  // radius(km),mass(solarmasses),energy density(MeV/fm^3),pressure(MeV/fm^3)
    fprintf(fp,"%le\t%le\t%le\t%le\t%le\n"
      ,epsi_zero* BAG, r/ 1.0e18, mass_solar, rTOV/ 1.0e18, massTOV_solar);
    
    epsi_zero= epsi_zero+ epsi_step;   
  }

  printf("\nProgram complete.\n\n");
  return EXIT_SUCCESS;
}

// FUNCTIONS ******************************************************************
//Epsilon differential equation using Newtonian physics
double dedr(double epsi, double mass, double r) {
  return ((-3.0* K* epsi* mass)/(pow(r,2.0)));
}

//mass differential equation
double dmdr(double epsi, double r) {
  return (4.0* M_PI* pow(r,2.0)* epsi);
}

//Epsilon differential equation using TOV formula
double dedrTOV(double epsi, double mass, double r) {
  return (-4.0* K* (epsi - BAG)*
  (mass + ((4.0* M_PI* pow(r,3.0))/ 3.0)* (epsi- 4.0* BAG))/
  (pow(r,2.0)- 2.0* K* mass* r));
}
