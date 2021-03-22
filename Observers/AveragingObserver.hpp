//
// Created by Nikita Kruk on 23.07.20.
//

#ifndef SPSSLANGEVININTEGRATION_AVERAGINGOBSERVER_HPP
#define SPSSLANGEVININTEGRATION_AVERAGINGOBSERVER_HPP

#include "../Definitions.hpp"

#include <vector>

class AveragingObserver
{
 public:

  AveragingObserver(int number_of_particles, int number_of_state_variables);
  ~AveragingObserver();

  void operator()(const std::vector<Real> &system_state, Real t);
  void SaveAveragedSystemState(Real v_0,
                               Real xi,
                               Real sigma,
                               Real rho,
                               Real alpha,
                               Real D_phi,
                               int trial);

 private:

  int n_;
  int s_;
  std::vector<Real> accumulated_system_state_;
  int number_of_accumulated_steps_;

};

#endif //SPSSLANGEVININTEGRATION_AVERAGINGOBSERVER_HPP
