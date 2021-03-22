//
// Created by Nikita Kruk on 14.02.20.
//

#ifndef SPPKURAMOTOWITHINERTIAODEINTEGRATION_THREADSHAREDMEMORY_HPP
#define SPPKURAMOTOWITHINERTIAODEINTEGRATION_THREADSHAREDMEMORY_HPP

#include "../Definitions.hpp"
#include "Thread.hpp"

#include <mpi.h>
#include <vector>
#include <unordered_map>

class ThreadSharedMemory : public Thread
{
 public:

  ThreadSharedMemory(int argc, char **argv, int number_of_particles, int number_of_state_variables);
  ~ThreadSharedMemory();

  int GetRootRank();
  int GetSharedRootRank();
  MPI_Comm &GetSharedCommunicator();
  MPI_Win &GetWindow(const std::string &window_name);

  bool IsRoot() const override;
  bool IsSharedRoot();

  void SynchronizeVector(std::vector<Real> &vec) override;
  void SynchronizeVector(Real *const vec, long size);
  void SynchronizeVectorThoughBuffer(std::vector<Real> &vec, std::vector<Real> &buf) override;
  void SynchronizeVectorThroughoutClusters(Real *const vec);
  void SynchronizeVectorThroughoutClusters(int *const vec);
  void BroadcastVector(std::vector<Real> &vec) override;
  void BroadcastVector(Real *const vec, long size);
  void SynchronizePrelinkedList(std::vector<int> &pre_linked_list) const override;

  void AllocateSharedWindow(MPI_Aint size, Real *&vec, const std::string &window_name);
  void AllocateSharedWindow(MPI_Aint size, int *&vec, const std::string &window_name);
  void FreeSharedWindow(const std::string &window_name);

 private:

  int shared_rank_;
  std::vector<int> shared_ranks_;
  int shared_root_rank_;
  std::vector<int> shared_root_rank_indexes_;
  int number_of_shared_mpich_threads_; // in a shared memory communicator
  std::vector<int> numbers_of_shared_mpich_threads_; // in each shared memory communicator
  MPI_Comm shared_communicator_;
  std::unordered_map<std::string, MPI_Win> windows_;
};

#endif //SPPKURAMOTOWITHINERTIAODEINTEGRATION_THREADSHAREDMEMORY_HPP
