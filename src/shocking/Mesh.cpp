#include <iostream>     /* cerr */
#include <stdlib.h>     /* atoi */

#include "Mesh.h"


Mesh::Mesh(const std::vector< std::string > l_input_file)
        : m_Xinitial(0.0),
          m_Xfinal(0.0),
          m_dX(0.0),
          m_dX_print(0.0)
{
    // Read mesh from input file
    inputFileParse(l_input_file);

    // Print mesh info
    printMeshInfo();
}

// --------------------------------------------------------

double Mesh::getXinitial(){return m_Xinitial; }
double Mesh::getXfinal(){return m_Xfinal;}
double Mesh::getdX(){return m_dX;}
double Mesh::getdXprint(){return m_dX_print;}

// --------------------------------------------------------

void Mesh::printMeshInfo()
{
  // Prints info about the mesh
  std::cout << "Mesh parameters:" << std::endl;
  std::cout << "    X_initial [m]: " << m_Xinitial << std::endl;
  std::cout << "    X_final [m]:   " << m_Xfinal   << std::endl;
  std::cout << "    dX [m]:        " << m_dX       << std::endl;
  std::cout << std::endl;
}

// --------------------------------------------------------

void Mesh::inputFileParse(const std::vector< std::string > l_input_file)
{
  // Parses the input file lines

  for(size_t id_l = 0; id_l < l_input_file.size(); ++id_l) {

    // X_initial
    if( l_input_file.at(id_l).compare("Mesh X_init:") == 0) {
      m_Xinitial = atof(l_input_file.at(id_l+1).c_str());
    }

    // X_final
    if( l_input_file.at(id_l).compare("Mesh X_final:") == 0) {
      m_Xfinal = atof(l_input_file.at(id_l+1).c_str());
    }

    // dX
    if( l_input_file.at(id_l).compare("Mesh dX:") == 0) {
      m_dX = atof(l_input_file.at(id_l+1).c_str());
    }
  }

}

// --------------------------------------------------------

void Mesh::errorMeshNotSetProperly()
{
  // Check if the mesh was set properly
  if(m_Xfinal <= m_Xinitial){
    std::cerr << " ATTENTION: X_final must be bigger than X_initial!"
              << " Check the input file!" << std::endl;
    exit(1);
  }

  if(m_dX <= 0){
    std::cerr << " ATTENTION: dX must be positive nonzero!"
              << " Check the input file!" << std::endl;
    exit(1);
  } 
}
