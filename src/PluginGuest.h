#ifndef PLUGINGUEST_H
#define PLUGINGUEST_H

#include "mpi.h"
#include <string>
#include <vector>

/** @example plugin-example-1.cpp */
/*!
 * \brief Supports the communication of a plugin that has been launched via MPI_Comm_spawn.
 * The communication is driven by the host, and follows the following convention.
 *   -# All communication is via printable character strings sent as MPI_CHAR using MPI_Send / MPI_Recv
 *   -# Strings are sent from host to guest by first sending the length of the string as an MPI_INT, then
 * 	the string itself.
 *   -# Strings are sent from guest to host by sending the length, then listening for the return
 * 	of an MPI_INT answer. If non-zero, the string is sent.
 *   -# The host sends an initial string containing the name of the host program, its
 * 	version number, and the desired file into which the guest's standard output should be redirected.
 * 	It then waits for strings to be sent to it, which, on receipt, it interprets and
 * 	obeys according to further conventions established between the two programs.
 *   -# The guest indicates termination by sending an empty string.
 *
 *   Example of use: @include plugin-example-1.cpp
 */
class PluginGuest {
public:
  /*!
   * \brief Construct a plugin guest instance
   * \param host is the name that the host program is expected to offer to identify itself. If the parameter is not given, then no check is performed.
   * \param world is the world MPI communicator. This function is collective across all processes in this communicator.
   */
  PluginGuest(std::string host="",MPI_Comm world=MPI_COMM_WORLD);
  ~PluginGuest() { close(); }
  bool active() const { return m_active;} ///< Whether the plugin is active or not
  std::string receive() const; ///< Receive a string from the host.
  /*!
   * \brief Send a string to the host
   * \param value The string to be sent
   * \return Whether the string was accepted by the host
   */
  bool send(const std::string value) const;
  /*!
   * \brief Close down the plugin. Called by the class destructor.
   */
  void close();
private:
  MPI_Comm m_intercomm; ///< the intercommunicator to the parent
  int m_active; ///< whether this Molpro plugin framework is active
  int m_rank; ///< rank of process. All that matters is that one of the processes has m_rank=0
  std::string m_host; ///< the name of the host program
  std::string m_hostVersion; ///< the version of the host program
  std::string m_outputFile; ///< the file on which standard output is written
};

// pure C interface
extern "C" {
  void PluginGuestOpen(char* host);
  int PluginGuestActive();
  int PluginGuestMaster();
  int PluginGuestSend(const char* value);
  const char* PluginGuestReceive();
  void PluginGuestClose();
}


#endif // PLUGINGUEST_H
