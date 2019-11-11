#include "PluginGuest.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <stdexcept>
#include <sstream>
PluginGuest::PluginGuest(const std::string host, const MPI_Comm world)
  : m_world(world)
{
  MPI_Comm_rank(m_world,&m_rank);
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
          MPI_Comm_size(m_world,&size);
          std::cout << "Plugin for "<<m_host<<" version "<<m_hostVersion<<" running with "<<size<<" processes"<<std::endl;
        }
  }
  {
    MPI_Bcast(&m_active,1,MPI_INT,0,m_world);
    int length=m_host.size();
    MPI_Bcast(&length,1,MPI_INT,0,m_world);
    if (m_rank!=0) m_host.resize(length);
    MPI_Bcast(&m_host[0],length,MPI_CHAR,0,m_world);
  }
  {
    int length=m_hostVersion.size();
    MPI_Bcast(&length,1,MPI_INT,0,m_world);
    if (m_rank!=0) m_hostVersion.resize(length);
    MPI_Bcast(&m_hostVersion[0],length,MPI_CHAR,0,m_world);
  }
}

bool PluginGuest::send(const std::string value) const
{
  if (! m_active) return true;
  int answer;
  int length=value.size();
  MPI_Bcast(&length,1,MPI_INT,0,m_world);
  if (m_rank==0)
    MPI_Send(&length,1,MPI_INT,0,0,m_intercomm);
  if (!length) return true;
  if (m_rank==0) {
    MPI_Status status;
    MPI_Recv(&answer,1,MPI_INT,0,0,m_intercomm,&status); // receive accept or reject
    if (answer) // 'yes' answer received
      MPI_Send(value.c_str(),length,MPI_CHAR,0,1,m_intercomm);
  }
  MPI_Bcast(&answer,1,MPI_INT,0,m_world);
  return answer>0 || length==0 ;
}

std::string PluginGuest::receive() const
{
  if (! m_active) return "";
  MPI_Status status;
  int length;
  if (m_rank==0)
    MPI_Recv(&length,1,MPI_INT,0,0,m_intercomm,&status);
  MPI_Bcast(&length,1,MPI_INT,0,m_world);
  if (length==0) throw std::logic_error("plugin request has failed");
  std::string result(length,' ');
  if (m_rank==0)
    MPI_Recv(&result[0],length,MPI_CHAR,0,1,m_intercomm,&status);
  MPI_Bcast(&result[0],length,MPI_CHAR,0,m_world);
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
// 24.01.2019 K.G. I had to change the constructor here to the default one (which does
// the same as guest=nullptr) because some older implementations do not support the
// construction from nullptr
static std::shared_ptr<PluginGuest> guest{};
void PluginGuestOpen(const char* host) { guest = std::make_shared<PluginGuest>(std::string(host)); }
int PluginGuestActive() { return guest!=nullptr && guest->active() ? 1 : 0;}
int PluginGuestSend(const char* value) { return guest->send(std::string(value)) ? 1 : 0 ; }
void PluginGuestReceive(char * string, size_t length) { 
  std::string recstring = guest->receive();
  if ( length <= recstring.size() ) throw std::length_error(recstring+" path is too long");
  strcpy(string,recstring.c_str());}
void PluginGuestClose() { guest->close(); }
