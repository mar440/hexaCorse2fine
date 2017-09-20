#include <stdio.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkHexahedron.h>
#include <vtkCellArray.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vector>
#include <math.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkGenericCell.h>
#include <vtkPolyData.h>
#include <vtkIdList.h>
#include <algorithm>
#include <string>


// To extract a sub part
#include <vtkThreshold.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelection.h>
#include <vtkInformation.h>
#include <vtkUnstructuredGridGeometryFilter.h>
#include <vtkAppendFilter.h>
#include <set>
#include <stdexcept>
#include <vtkConnectivityFilter.h>
//#include <vtkIdTypeArray.h>
#include <iostream>

#define VTK_CREATE(type,name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

//
using namespace std;


void printVecVec(vector < vector < int> > &v){

    for (int i = 0; i < v.size();i++){
        for( int j = 0; j < v[i].size();j++){
            cout << v[i][j] << " ";
        }
        cout << endl;
    }
}


class Line
{
    public:
        int a;  // first point
        int b;  // second point
        int i;  // order number
        int s;
        bool reverted;
        void setVal(int,int,int);

        Line();
        Line(int,int,int);
        ~Line();

};

void Line::setVal(int i_, int a_, int b_) {
    this->i = i_;
    this->s = 0;

    if (a_ > b_){
        this->a = b_;
        this->b = a_;
        this->reverted = true;
    }
    else{
        this->a = a_;
        this->b = b_;
        this->reverted = false;
    }
}

Line::Line(){
}

Line::Line(int i_, int a_, int b_){
    Line::setVal(i_,a_,b_);
}

Line::~Line(){
}

// -------------------------------------------
class Face
{
    public:
        int a;
        int b;
        int c;
        int d;
        int i;
        int s;
        bool reverted;
        void setVal(int,int,int,int,int);

        Face();
        Face(int,int,int,int,int);
        ~Face();

};

void Face::setVal(int i_, int a_, int b_, int c_, int d_) {
    this->i = i_;
    this->s = 0;

    vector <int> v_({a_,b_,c_,d_});
//    sort(v_.begin(),v_.end());

    this->a = v_[0];
    this->b = v_[1];
    this->c = v_[2];
    this->d = v_[3];

}

Face::Face(){
}

Face::Face(int i_, int a_, int b_, int c_, int d_){
    Face::setVal(i_,a_,b_,c_,d_);
}

Face::~Face(){
}


bool compare_a(Line x, Line y){
    return x.a < y.a;
}
bool compare_b(Line x, Line y){
    return x.b < y.b;
}


bool Fcompare_a(Face x, Face y){
    return x.a < y.a;
}

bool Fcompare_b(Face x, Face y){
    return x.b < y.b;
}


bool Fcompare_c(Face x, Face y){
    return x.c < y.c;
}

bool Fcompare_d(Face x, Face y){
    return x.d < y.d;
}

void printLines(vector <Line> x){
    for (int i = 0; i < x.size(); i++){
        cout << "order: " << x[i].i <<",  ";
        cout <<  x[i].a <<" "<< x[i].b <<", ";
        cout <<  "s   " << x[i].s << endl;
    }
}

void printFaces(vector <Face> x){
    for (int i = 0; i < x.size(); i++){
        cout << "order: " << x[i].i <<",  ";
        cout <<  x[i].a <<", "<< x[i].b <<", "<< x[i].c <<", "<< x[i].d <<", ";
        cout <<  "s   " << x[i].s << endl;
    }
}

void sortLines(vector < Line > &v, int &nP){
    sort(v.begin(),v.end(),compare_a);

    int startInd = 0, endInd = 1;
    int tmpVecIprev = v[0].a;
    int nnz = v.size();
    for (int i = 1 ; i < nnz; i ++){
        if (v[i].a == tmpVecIprev){
            endInd++;
        }
        if (v[i].a != tmpVecIprev || (v[i].a == tmpVecIprev &&
                 i == nnz - 1))
        {
            sort(v.begin() + startInd,
                              v.begin() + (endInd  ),compare_b);
            startInd = i;
            endInd = i + 1;
        }
        tmpVecIprev = v[i].a;
    }


    v[0].s = nP;
    for (int i = 1; i < nnz; i++){
        if (v[i-1].a == v[i].a && v[i-1].b == v[i].b){
            // nothing
        }
        else{
            nP++;
        }
        v[i].s = nP;
    }
    nP++;
}


void sortFaces(vector < Face > &v, int &nP){
    sort(v.begin(),v.end(),Fcompare_a);

    int startInd = 0, endInd = 1;
    int tmpVecIprev = v[0].a;
    int nnz = v.size();
    for (int i = 1 ; i < nnz; i ++){
        if (v[i].a == tmpVecIprev){
            endInd++;
        }
        if (v[i].a != tmpVecIprev || (v[i].a == tmpVecIprev &&
                 i == nnz - 1))
        {
            sort(v.begin() + startInd,
                              v.begin() + (endInd  ),Fcompare_b);
            startInd = i;
            endInd = i + 1;
        }
        tmpVecIprev = v[i].a;
    }


    startInd = 0, endInd = 1;
    tmpVecIprev = v[0].b;

    for (int i = 1 ; i < nnz; i ++){
        if (v[i].b == tmpVecIprev){
            endInd++;
        }
        if (v[i].b != tmpVecIprev || (v[i].b == tmpVecIprev &&
                 i == nnz - 1))
        {
            sort(v.begin() + startInd,
                              v.begin() + (endInd  ),Fcompare_c);
            startInd = i;
            endInd = i + 1;
        }
        tmpVecIprev = v[i].b;
    }

    startInd = 0, endInd = 1;
    tmpVecIprev = v[0].c;

    for (int i = 1 ; i < nnz; i ++){
        if (v[i].c == tmpVecIprev){
            endInd++;
        }
        if (v[i].c != tmpVecIprev || (v[i].c == tmpVecIprev &&
                 i == nnz - 1))
        {
            sort(v.begin() + startInd,
                              v.begin() + (endInd  ),Fcompare_d);
            startInd = i;
            endInd = i + 1;
        }
        tmpVecIprev = v[i].c;
    }


    v[0].s = nP;
    for (int i = 1; i < nnz; i++){
        if (v[i-1].a == v[i].a && v[i-1].b == v[i].b &&
                v[i-1].c == v[i].c && v[i-1].d == v[i].d){
            // nothing
        }
        else{
            nP++;
        }
        v[i].s = nP;
    }



}





int main(int argc, char *argv[])
{

    string filename = "../output.vtu";


    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    vtkSmartPointer<vtkGenericCell> faceCell = vtkSmartPointer<vtkGenericCell>::New();
    //read all the data from the file
    vtkSmartPointer<vtkXMLUnstructuredGridReader> mesh=
            vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    mesh->SetFileName(filename.c_str());
    mesh->Update();

    const int nCellsOrig  = mesh->GetNumberOfCells();
    int nCells  = nCellsOrig;

    const int nPointsOrig = mesh->GetNumberOfPoints();
    int nPoints = nPointsOrig;

    cout << "nCellsOrig  = " << nCellsOrig << endl;
    cout << "nPointsOrig = " << nPointsOrig << endl;

    // only cube considered
    vector < vector < int> > edges;
    vector < vector < int> > faces;

    for (int i = 0 ; i < nCellsOrig; i++){

        mesh->GetOutput()->GetCell(i,cell);

        int i0 = cell->GetPointId(0);
        int i1 = cell->GetPointId(1);
        int i2 = cell->GetPointId(2);
        int i3 = cell->GetPointId(3);
        int i4 = cell->GetPointId(4);
        int i5 = cell->GetPointId(5);
        int i6 = cell->GetPointId(6);
        int i7 = cell->GetPointId(7);

//          7-------6
//       4------5   :
//       :  :   :   :
//       :  :   :   :
//       :  3---:---2
//       0------1

        vector <int> v_({i0,i1,i5,i4});
        sort(v_.begin(),v_.end());
        faces.push_back(v_);

        v_ = {i1,i2,i6,i5};
        sort(v_.begin(),v_.end());
        faces.push_back(v_);

        v_ = {i2,i3,i7,i6};
        sort(v_.begin(),v_.end());
        faces.push_back(v_);

        v_ = {i3,i0,i4,i7};
        sort(v_.begin(),v_.end());
        faces.push_back(v_);

        v_ = {i4,i5,i6,i7};
        sort(v_.begin(),v_.end());
        faces.push_back(v_);

        v_= {i3,i2,i1,i0};
        sort(v_.begin(),v_.end());
        faces.push_back(v_);


        edges.push_back({i0,i1});
        edges.push_back({i1,i2});
        edges.push_back({i2,i3});
        edges.push_back({i3,i0});

        edges.push_back({i4,i5});
        edges.push_back({i5,i6});
        edges.push_back({i6,i7});
        edges.push_back({i7,i4});

        edges.push_back({i0,i4});
        edges.push_back({i1,i5});
        edges.push_back({i2,i6});
        edges.push_back({i3,i7});

    }





    vector <  Line  > edg;
    edg.resize(edges.size());
    for (int i = 0 ; i < edges.size() ; i++){
        edg[i].setVal(i,edges[i][0], edges[i][1]);
    }
    vector <  Line  > tmp_edg= edg;
    int nPCurrent = nPointsOrig;
    sortLines(tmp_edg,nPCurrent);
    //printLines(tmp_edg);


    for (int i = 0; i < edges.size(); i++){
        edges[tmp_edg[i].i].push_back(tmp_edg[i].s);
    }

    //printVecVec(edges);




    vector <  Face > fc;
    fc.resize(faces.size());
    for (int i = 0 ; i < fc.size() ; i++){
        fc[i].setVal(i,faces[i][0], faces[i][1], faces[i][2], faces[i][3]);
    }
    vector <  Face > tmp_fc= fc;
    sortFaces(tmp_fc,nPCurrent);
    //printFaces(tmp_fc);


    for (int i = 0; i < faces.size(); i++){
        faces[tmp_fc[i].i].push_back(tmp_fc[i].s);
    }


    //printVecVec(faces);


    vtkSmartPointer<vtkPoints> newPoints = vtkSmartPointer<vtkPoints>::New();
    int estim = nCellsOrig * 27;
    newPoints->SetNumberOfPoints(estim);

    int iNew[27];
    double ijPoint0[3], ijPoint1[3];
    double x_,y_,z_;


    vtkSmartPointer<vtkHexahedron> hexa = vtkSmartPointer<vtkHexahedron>::New();
// Cell array - connectivity
    vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();


    for (int i = 0; i < nCellsOrig ; i++)
    {
        mesh->GetOutput()->GetCell(i,cell);

//      firt 8 nodes are original nodes of parent hexa8
        for (int j = 0; j < 8 ; j++){
            iNew[j] = cell->GetPointId(j);
            mesh->GetOutput()->GetPoint(iNew[j],ijPoint0);
            newPoints->SetPoint(iNew[j],ijPoint0[0],ijPoint0[1],ijPoint0[2]);
        }
//      next 12 nodes come from splinting all edges
        iNew[ 8] = edges[12 * i +  0][2];
        mesh->GetOutput()->GetPoint(edges[12 * i +  0][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i +  0][1],ijPoint1);
        newPoints->SetPoint(iNew[8], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                        0.5 * (ijPoint0[1] + ijPoint1[1]),
                                        0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[ 9] = edges[12 * i +  1][2];
        mesh->GetOutput()->GetPoint(edges[12 * i +  1][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i +  1][1],ijPoint1);
        newPoints->SetPoint(iNew[9], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                        0.5 * (ijPoint0[1] + ijPoint1[1]),
                                        0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[10] = edges[12 * i +  2][2];
        mesh->GetOutput()->GetPoint(edges[12 * i +  2][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i +  2][1],ijPoint1);
        newPoints->SetPoint(iNew[10], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                         0.5 * (ijPoint0[1] + ijPoint1[1]),
                                         0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[11] = edges[12 * i +  3][2];
        mesh->GetOutput()->GetPoint(edges[12 * i +  3][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i +  3][1],ijPoint1);
        newPoints->SetPoint(iNew[11], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                         0.5 * (ijPoint0[1] + ijPoint1[1]),
                                         0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[12] = edges[12 * i +  4][2];
        mesh->GetOutput()->GetPoint(edges[12 * i +  4][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i +  4][1],ijPoint1);
        newPoints->SetPoint(iNew[12], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                         0.5 * (ijPoint0[1] + ijPoint1[1]),
                                         0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[13] = edges[12 * i +  5][2];
        mesh->GetOutput()->GetPoint(edges[12 * i +  5][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i +  5][1],ijPoint1);
        newPoints->SetPoint(iNew[13], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                         0.5 * (ijPoint0[1] + ijPoint1[1]),
                                         0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[14] = edges[12 * i +  6][2];
        mesh->GetOutput()->GetPoint(edges[12 * i +  6][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i +  6][1],ijPoint1);
        newPoints->SetPoint(iNew[14], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                         0.5 * (ijPoint0[1] + ijPoint1[1]),
                                         0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[15] = edges[12 * i +  7][2];
        mesh->GetOutput()->GetPoint(edges[12 * i +  7][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i +  7][1],ijPoint1);
        newPoints->SetPoint(iNew[15], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                         0.5 * (ijPoint0[1] + ijPoint1[1]),
                                         0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[16] = edges[12 * i +  8][2];
        mesh->GetOutput()->GetPoint(edges[12 * i +  8][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i +  8][1],ijPoint1);
        newPoints->SetPoint(iNew[16], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                         0.5 * (ijPoint0[1] + ijPoint1[1]),
                                         0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[17] = edges[12 * i +  9][2];
        mesh->GetOutput()->GetPoint(edges[12 * i +  9][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i +  9][1],ijPoint1);
        newPoints->SetPoint(iNew[17], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                         0.5 * (ijPoint0[1] + ijPoint1[1]),
                                         0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[18] = edges[12 * i + 10][2];
        mesh->GetOutput()->GetPoint(edges[12 * i + 10][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i + 10][1],ijPoint1);
        newPoints->SetPoint(iNew[18], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                         0.5 * (ijPoint0[1] + ijPoint1[1]),
                                         0.5 * (ijPoint0[2] + ijPoint1[2]));

        iNew[19] = edges[12 * i + 11][2];
        mesh->GetOutput()->GetPoint(edges[12 * i + 11][0],ijPoint0);
        mesh->GetOutput()->GetPoint(edges[12 * i + 11][1],ijPoint1);
        newPoints->SetPoint(iNew[19], 0.5 * (ijPoint0[0] + ijPoint1[0]),
                                         0.5 * (ijPoint0[1] + ijPoint1[1]),
                                         0.5 * (ijPoint0[2] + ijPoint1[2]));

//      next 6 nodes come from each face (center)
        iNew[20] = faces[6 * i +  0][4];
        mesh->GetOutput()->GetPoint(faces[6 * i +  0][0],ijPoint1);
        ijPoint0[0] = ijPoint1[0];
        ijPoint0[1] = ijPoint1[1];
        ijPoint0[2] = ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  0][1],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  0][2],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  0][3],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        newPoints->SetPoint(iNew[20], 0.25 * ijPoint0[0],
                                      0.25 * ijPoint0[1],
                                      0.25 * ijPoint0[2]);


        iNew[21] = faces[6 * i +  1][4];
        mesh->GetOutput()->GetPoint(faces[6 * i +  1][0],ijPoint1);
        ijPoint0[0] = ijPoint1[0];
        ijPoint0[1] = ijPoint1[1];
        ijPoint0[2] = ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  1][1],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  1][2],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  1][3],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        newPoints->SetPoint(iNew[21], 0.25 * ijPoint0[0],
                                         0.25 * ijPoint0[1],
                                         0.25 * ijPoint0[2]);


        iNew[22] = faces[6 * i +  2][4];
        mesh->GetOutput()->GetPoint(faces[6 * i +  2][0],ijPoint1);
        ijPoint0[0] = ijPoint1[0];
        ijPoint0[1] = ijPoint1[1];
        ijPoint0[2] = ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  2][1],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  2][2],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  2][3],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        newPoints->SetPoint(iNew[22], 0.25 * ijPoint0[0],
                                         0.25 * ijPoint0[1],
                                         0.25 * ijPoint0[2]);


        iNew[23] = faces[6 * i +  3][4];
        mesh->GetOutput()->GetPoint(faces[6 * i +  3][0],ijPoint1);
        ijPoint0[0] = ijPoint1[0];
        ijPoint0[1] = ijPoint1[1];
        ijPoint0[2] = ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  3][1],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  3][2],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  3][3],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        newPoints->SetPoint(iNew[23], 0.25 * ijPoint0[0],
                                         0.25 * ijPoint0[1],
                                         0.25 * ijPoint0[2]);


        iNew[24] = faces[6 * i +  4][4];
        mesh->GetOutput()->GetPoint(faces[6 * i +  4][0],ijPoint1);
        ijPoint0[0] = ijPoint1[0];
        ijPoint0[1] = ijPoint1[1];
        ijPoint0[2] = ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  4][1],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  4][2],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  4][3],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        newPoints->SetPoint(iNew[24], 0.25 * ijPoint0[0],
                                         0.25 * ijPoint0[1],
                                         0.25 * ijPoint0[2]);


        iNew[25] = faces[6 * i +  5][4];
        mesh->GetOutput()->GetPoint(faces[6 * i +  5][0],ijPoint1);
        ijPoint0[0] = ijPoint1[0];
        ijPoint0[1] = ijPoint1[1];
        ijPoint0[2] = ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  5][1],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  5][2],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        mesh->GetOutput()->GetPoint(faces[6 * i +  5][3],ijPoint1);
        ijPoint0[0] += ijPoint1[0];
        ijPoint0[1] += ijPoint1[1];
        ijPoint0[2] += ijPoint1[2];
        newPoints->SetPoint(iNew[25], 0.25 * ijPoint0[0],
                                         0.25 * ijPoint0[1],
                                         0.25 * ijPoint0[2]);


//      next 1 node comes from center of parent element
        iNew[26] = nPCurrent;

        ijPoint0[0] = 0; ijPoint0[1] = 0; ijPoint0[2] = 0;
        for (int j = 0; j < 8; j++ ){
            mesh->GetOutput()->GetPoint(iNew[j],ijPoint1);
            ijPoint0[0] += ijPoint1[0];
            ijPoint0[1] += ijPoint1[1];
            ijPoint0[2] += ijPoint1[2];
        }
        newPoints->SetPoint(iNew[26], 0.125 * ijPoint0[0],
                                         0.125 * ijPoint0[1],
                                         0.125 * ijPoint0[2]);
        nPCurrent++;


        int ind00[8][8] = { { 0, 8,25,11,16,20,26,23},
                            { 8, 1, 9,25,20,17,21,26},
                            {11,25,10, 3,23,26,22,19},
                            {25, 9, 2,10,26,21,18,22},
                            {16,20,26,23, 4,12,24,15},
                            {20,17,21,26,12, 5,13,24},
                            {23,26,22,19,15,24,14, 7},
                            {26,21,18,22,24,13, 6,14}};

        int iii;
        for (int j = 0; j < 8; j++)
 //       int j = 1;
        {
            for (int k = 0; k < 8; k++){
                iii = iNew[ind00[j][k]];
                hexa->GetPointIds()->SetId(k, iii);
            }
            cellArray->InsertNextCell(hexa);
        }


    }


    newPoints->SetNumberOfPoints(nPCurrent);

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
        vtkSmartPointer<vtkUnstructuredGrid>::New();
      unstructuredGrid->SetPoints(newPoints);
      unstructuredGrid->SetCells(VTK_HEXAHEDRON, cellArray);

  vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
    vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  writer->SetFileName("../output.vtu");

  bool asciiOrBinaryVtu = false;
  if (asciiOrBinaryVtu){
    writer->SetDataModeToAscii();
  }
  else{
    writer->SetDataModeToBinary();
  }


#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(unstructuredGrid);
#else
  writer->SetInputData(unstructuredGrid);
#endif
  writer->Write();

    cout << " edg size          :   " << edg.size() << endl;
    cout << " nCellsOrig * 6 * 4:   " << nCellsOrig * 12 << endl;
    cout << " nCellsOrig        :   " << nCellsOrig  << endl;

    return 0;

}
