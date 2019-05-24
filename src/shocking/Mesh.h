#ifndef MESH_H
#define MESH_H

#include <string>
#include <vector>

class Mesh {
public:
    Mesh(const std::vector< std::string > l_input_file);
    ~Mesh(){}

    double getXinitial();
    double getXfinal();
    double getdX();
    double getdXprint();

private:
    double m_Xinitial;
    double m_Xfinal;
    double m_dX;
    double m_dX_print;

    void inputFileParse(const std::vector< std::string > l_input_file);
    void errorMeshNotSetProperly();
    void printMeshInfo();

};

#endif /* MESH_H */
