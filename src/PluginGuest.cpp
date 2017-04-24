#include "PluginGuest.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <stdexcept>
#include <sstream>
PluginGuest::PluginGuest(std::string host, MPI_Comm world)
{
  MPI_Comm_rank(world,&m_rank);
  m_active=0;
  MPI_Comm_get_parent(&m_intercomm);
  if (m_rank==0 && m_intercomm != MPI_COMM_NULL) {
      int length;
      MPI_Status status;
      // expect plugin server to identify itself
      MPI_Recv(&length,1,MPI_INT,0,0,m_intercomm,&status);
      std::string id(' ',length);
      MPI_Recv(&id[0],length,MPI_CHAR,0,1,m_intercomm,&status);
      std::stringstream ss(id);
      ss >> m_host;
      m_active = host=="" || host == m_host;
      if (m_active) {
          ss >> m_hostVersion;
          std::string m_outputFile;
          ss >> m_outputFile;
          if (m_outputFile != "") {
              fflush(stdout);
              if (freopen(&m_outputFile[0],"w",stdout)==NULL)
                perror("Error in opening plugin output file");
            }
          int size;
          MPI_Comm_size(world,&size);
          std::cout << "Plugin for "<<m_host<<" version "<<m_hostVersion<<" running with "<<size<<" processes"<<std::endl;
        }
    }
  MPI_Bcast(&m_active,1,MPI_INT,0,world);
}

bool PluginGuest::send(const std::string value) const
{
  if (m_rank>0 || ! m_active) return true;
  int length=value.size();
  MPI_Send(&length,1,MPI_INT,0,0,m_intercomm);
  if (!length) return true;
  int answer;
  MPI_Status status;
  MPI_Recv(&answer,1,MPI_INT,0,0,m_intercomm,&status); // receive accept or reject
  if (answer) // 'yes' answer received
    MPI_Send(value.c_str(),length,MPI_CHAR,0,1,m_intercomm);
  return answer>0;
}

std::string PluginGuest::receive() const
{
  if (m_rank>0 || ! m_active) return "";
  MPI_Status status;
  int length;
  MPI_Recv(&length,1,MPI_INT,0,0,m_intercomm,&status);
  if (length==0) throw std::logic_error("plugin request has failed");
  std::string result(length,' ');
  MPI_Recv(&result[0],length,MPI_CHAR,0,1,m_intercomm,&status);
  return result;
}

void PluginGuest::close()
{
  fflush(stdout);
  send("");
  m_active=false;
}

// pure C interface
#include <memory>
static std::shared_ptr<PluginGuest> guest = nullptr;
void PluginGuestOpen(char* host) { guest = std::make_shared<PluginGuest>(std::string(host)); }
int PluginGuestActive() { return guest!=nullptr && guest->active() ? 1 : 0;}
int PluginGuestMaster() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  return rank==0;
}
int PluginGuestSend(const char* value) { return guest->send(std::string(value)) ? 1 : 0 ; }
const char* PluginGuestReceive() { return guest->receive().c_str(); }
void PluginGuestClose() { guest->close(); }
