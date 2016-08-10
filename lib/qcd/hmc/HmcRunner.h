    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/qcd/hmc/HmcRunner.h

    Copyright (C) 2015

Author: paboyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef HMC_RUNNER
#define HMC_RUNNER

namespace Grid{
  namespace QCD{


template<class Gimpl>
class NerscHmcRunnerTemplate {
public:

  INHERIT_GIMPL_TYPES(Gimpl);

  enum StartType_t { ColdStart, HotStart, TepidStart, CheckpointStart };

  ActionSet<GaugeField> TheAction;

  GridCartesian         * UGrid   ;
  GridCartesian         * FGrid   ;
  GridRedBlackCartesian * UrbGrid ;
  GridRedBlackCartesian * FrbGrid ;

  virtual void BuildTheAction (int argc, char **argv) = 0; // necessary?

  
  void Run (int argc, char  **argv){

    StartType_t StartType = HotStart;

    std::string arg;

    if( GridCmdOptionExists(argv,argv+argc,"--StartType") ){
      arg = GridCmdOptionPayload(argv,argv+argc,"--StartType");
      if ( arg == "HotStart" ) { StartType = HotStart; }
      else if ( arg == "ColdStart" ) { StartType = ColdStart; }
      else if ( arg == "TepidStart" ) { StartType = TepidStart; }
      else if ( arg == "CheckpointStart" ) { StartType = CheckpointStart; }
      else assert(0);
    }

    int StartTraj = 0;
    if( GridCmdOptionExists(argv,argv+argc,"--StartTrajectory") ){
      arg= GridCmdOptionPayload(argv,argv+argc,"--StartTrajectory");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg,ivec);
      StartTraj = ivec[0];
    }    

    int NumTraj = 1;
    if( GridCmdOptionExists(argv,argv+argc,"--Trajectories") ){
      arg= GridCmdOptionPayload(argv,argv+argc,"--Trajectories");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg,ivec);
      NumTraj = ivec[0];
    }

    int NumThermalizations = 10;
    if( GridCmdOptionExists(argv,argv+argc,"--Thermalizations") ){
      arg= GridCmdOptionPayload(argv,argv+argc,"--Thermalizations");
      std::vector<int> ivec(0);
      GridCmdOptionIntVector(arg,ivec);
      NumThermalizations = ivec[0];
    }


    GridSerialRNG    sRNG;
    GridParallelRNG  pRNG(UGrid);
    GaugeField  U(UGrid); // change this to an extended field (smearing class)

    std::vector<int> SerSeed({1,2,3,4,5});
    std::vector<int> ParSeed({6,7,8,9,10});

    
    // Create integrator, including the smearing policy
    // Smearing policy
#if 0
    std::cout << GridLogDebug << " Creating the Stout class\n";
    double rho = 0.1; // smearing parameter, now hardcoded
    int Nsmear = 1;   // number of smearing levels
    Smear_Stout<Gimpl> Stout(rho);
    std::cout << GridLogDebug << " Creating the SmearedConfiguration class\n";
    SmearedConfiguration<Gimpl> SmearingPolicy(UGrid, Nsmear, Stout);
    std::cout << GridLogDebug << " done\n";
#else 
    NoSmearing<Gimpl> SmearingPolicy;
#endif 
    //////////////

    // change here to change the algorithm
    typedef MinimumNorm2<Gimpl, NoSmearing<Gimpl> >  IntegratorType;
    IntegratorParameters MDpar(20);
    IntegratorType MDynamics(UGrid, MDpar, TheAction, SmearingPolicy);

    
    // Checkpoint strategy
    //    NerscHmcCheckpointer<Gimpl> Checkpoint(std::string("ckpoint_lat"),std::string("ckpoint_rng"),1);
    PlaquetteLogger<Gimpl>      PlaqLog(std::string("plaq"));

    HMCparameters HMCpar;
    HMCpar.StartTrajectory   = StartTraj;
    HMCpar.Trajectories      = NumTraj;
    HMCpar.NoMetropolisUntil = NumThermalizations;
    

    if ( StartType == HotStart ) {
      // Hot start
      HMCpar.MetropolisTest = true;
      sRNG.SeedFixedIntegers(SerSeed);
      pRNG.SeedFixedIntegers(ParSeed);
      SU3::HotConfiguration(pRNG, U,Ndim);
    } else if ( StartType == ColdStart ) { 
      // Cold start
      HMCpar.MetropolisTest = true;
      sRNG.SeedFixedIntegers(SerSeed);
      pRNG.SeedFixedIntegers(ParSeed);
      SU3::ColdConfiguration(pRNG, U,Ndim);
    } else if ( StartType == TepidStart ) {       
      // Tepid start
      HMCpar.MetropolisTest = true;
      sRNG.SeedFixedIntegers(SerSeed);
      pRNG.SeedFixedIntegers(ParSeed);
      SU3::TepidConfiguration(pRNG, U,Ndim);
    } else if ( StartType == CheckpointStart ) { 
      HMCpar.MetropolisTest = true;
      // CheckpointRestart
      //      Checkpoint.CheckpointRestore(StartTraj, U, sRNG, pRNG);
    }

    // Attach the gauge field to the smearing Policy and create the fill the smeared set
    // notice that the unit configuration is singular in this procedure
    std::cout << GridLogMessage << "Filling the smeared set\n"; 
    SmearingPolicy.set_GaugeField(U);
    
    HybridMonteCarlo<GaugeField,IntegratorType>  HMC(HMCpar, MDynamics,sRNG,pRNG,U); 
    //    HMC.AddObservable(&Checkpoint);
    HMC.AddObservable(&PlaqLog);
    
    // Run it
    HMC.evolve();
    
  }
  
};

 typedef NerscHmcRunnerTemplate<PeriodicGimplR> NerscHmcRunner;
 typedef NerscHmcRunnerTemplate<PeriodicGimplF> NerscHmcRunnerF;
 typedef NerscHmcRunnerTemplate<PeriodicGimplD> NerscHmcRunnerD;

 typedef NerscHmcRunnerTemplate<PeriodicGimplR> PeriodicNerscHmcRunner;
 typedef NerscHmcRunnerTemplate<PeriodicGimplF> PeriodicNerscHmcRunnerF;
 typedef NerscHmcRunnerTemplate<PeriodicGimplD> PeriodicNerscHmcRunnerD;

 typedef NerscHmcRunnerTemplate<ConjugateGimplR> ConjugateNerscHmcRunner;
 typedef NerscHmcRunnerTemplate<ConjugateGimplF> ConjugateNerscHmcRunnerF;
 typedef NerscHmcRunnerTemplate<ConjugateGimplD> ConjugateNerscHmcRunnerD;

}}
#endif
