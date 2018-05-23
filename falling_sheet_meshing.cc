//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1307 $
//LIC//
//LIC// $LastChangedDate: 2018-01-18 11:30:14 +0000 (Thu, 18 Jan 2018) $
//LIC// 
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================

//Generic routines
#include "generic.h"

// Navier--Stokes equations
#include "navier_stokes.h"

// Poisson
#include "poisson.h"

// Get the mesh
#include "meshes/tetgen_mesh.h" 

// Local version of mesh
#include "refineable_tetgen_mesh.h"


// Do Navier Stokes (Fluid) or Poisson?
//#define DO_FLUID

using namespace std;

using namespace oomph;


//=============================================================
/// Namespace for problem parameters
//=============================================================
namespace Global_Parameters
{

 /// Reynolds number
 double Re = 0.0;

 /// (Half-)width of the box
 double Box_width = 1.5;

 /// (Half)height of the box
 double Box_length = 2.0;

}

//====================================================================
/// Disk in container
//====================================================================
template<class ELEMENT> 
class DiskInContainerProblem : public Problem
{

public:

 /// Constructor
 DiskInContainerProblem();
  
 /// Destructor (empty)
 ~DiskInContainerProblem()
  {
   //Delete the objects
   unsigned nh = Inner_boundary_pt.size();
   for(unsigned h=0;h<nh;++h)
    {
     delete Inner_boundary_pt[h];
    }
   delete Outer_boundary_pt;
  }

      
 /// Kill line visualiser
 void actions_before_adapt()
  {
   delete LV_pt;
   LV_pt=0;
  }

 /// Totally new mesh; build elements and apply boundary conditions
 void actions_after_adapt()
  {
   complete_problem_setup();
  }
 

 /// Build line visualiser
 void build_line_visualiser()
  {
   // How many points to you want?
   unsigned int  npt=1000;
   Vector<Vector<double> > coord_vec(npt);
   for (unsigned j=0;j<npt;j++)
    {
     coord_vec[j].resize(3);
     double coord=-Global_Parameters::Box_width+
      2.0*Global_Parameters::Box_width*
      double(j)/double(npt-1);
     coord_vec[j][0]=coord;
     coord_vec[j][1]=coord;
     coord_vec[j][2]=0.0;
    }

   // Setup line visualiser
   LV_pt=new LineVisualiser(Problem::mesh_pt(),
                            coord_vec);
   
  }



 /// Update the problem specs before solve: (Re)set boundary conditions
 void actions_before_newton_solve()
  { 
   // Sheet rises vertically
   unsigned ibound=Sheet_boundary_id;
   {
    // Loop over the nodes on boundary
    unsigned num_nod=mesh_pt()->nboundary_node(ibound);
    for (unsigned inod=0;inod<num_nod;inod++)
     {
      Node* nod_pt=mesh_pt()->boundary_node_pt(ibound,inod);
#ifdef DO_FLUID   
      nod_pt->set_value(2,1.0);
#else
      nod_pt->set_value(0,1.0);
#endif

     }
   }
  }

 /// Update the problem specs before solve (empty)
 void actions_after_newton_solve(){}

 
 //Access function for the specific mesh
 RefineableTetgenMesh<ELEMENT>* mesh_pt() 
  {
   return dynamic_cast<RefineableTetgenMesh<ELEMENT>*>(Problem::mesh_pt());
  }
 
 /// Doc the solution
 void doc_solution(const unsigned& nplot, DocInfo& doc_info);

private:
 
 /// Apply BCs and make elements functional
 void complete_problem_setup();

 /// Storage for the outer boundary object
 TetgenMeshFacetedSurface* Outer_boundary_pt;

 /// Inner boundary
 Vector<TetgenMeshFacetedSurface*> Inner_boundary_pt;

 /// Boundary ID: top
 unsigned Top_boundary_id;

 /// Boundary ID: bottom
 unsigned Bottom_boundary_id;

 /// Boundary ID: sheet
 unsigned Sheet_boundary_id;

 // The Line Visualiser.
 LineVisualiser* LV_pt;

};



//========================================================================
/// Constructor for DiskInContainer problem
//========================================================================
template<class ELEMENT>
DiskInContainerProblem<ELEMENT>::DiskInContainerProblem()
{ 

 //Add a steady time stepper
 this->add_time_stepper_pt(new Steady<0>);

 
 //Make the external box
 //---------------------
 const double box_width = Global_Parameters::Box_width;
 const double box_length = Global_Parameters::Box_length;
 Vector<Vector<double> > box_point(8);
 box_point[0].resize(3);
 box_point[0][0] = -box_width;
 box_point[0][1] = -box_width;
 box_point[0][2] = -box_length;

 box_point[1].resize(3);
 box_point[1][0] = -box_width;
 box_point[1][1] =  box_width;
 box_point[1][2] = -box_length;

 box_point[2].resize(3);
 box_point[2][0] = -box_width;
 box_point[2][1] =  box_width;
 box_point[2][2] =  box_length;

 box_point[3].resize(3);
 box_point[3][0] = -box_width;
 box_point[3][1] = -box_width;
 box_point[3][2] =  box_length;

 box_point[4].resize(3);
 box_point[4][0] =  box_width;
 box_point[4][1] = -box_width;
 box_point[4][2] = -box_length;

 box_point[5].resize(3);
 box_point[5][0] =  box_width;
 box_point[5][1] =  box_width;
 box_point[5][2] = -box_length;

 box_point[6].resize(3);
 box_point[6][0] =  box_width;
 box_point[6][1] =  box_width;
 box_point[6][2] =  box_length;

 box_point[7].resize(3);
 box_point[7][0] =  box_width;
 box_point[7][1] =  -box_width;
 box_point[7][2] =  box_length;

 Vector<Vector<unsigned> > box_facet(6);
 box_facet[0].resize(4);
 box_facet[0][0] = 0;
 box_facet[0][1] = 4;
 box_facet[0][2] = 7;
 box_facet[0][3] = 3;

 box_facet[1].resize(4);
 box_facet[1][0] = 4;
 box_facet[1][1] = 5;
 box_facet[1][2] = 6;
 box_facet[1][3] = 7;

 // top
 box_facet[2].resize(4);
 box_facet[2][0] = 3;
 box_facet[2][1] = 7;
 box_facet[2][2] = 6;
 box_facet[2][3] = 2;

 box_facet[3].resize(4);
 box_facet[3][0] = 6;
 box_facet[3][1] = 5;
 box_facet[3][2] = 1;
 box_facet[3][3] = 2;

 box_facet[4].resize(4);
 box_facet[4][0] = 3;
 box_facet[4][1] = 2;
 box_facet[4][2] = 1;
 box_facet[4][3] = 0;

 // bottom
 box_facet[5].resize(4);
 box_facet[5][0] = 5;
 box_facet[5][1] = 1;
 box_facet[5][2] = 0;
 box_facet[5][3] = 4;

 // Boundaries associated with facets
 Vector<unsigned> box_facet_boundary_id(6);
 box_facet_boundary_id[0] = 2;
 box_facet_boundary_id[1] = 3;
 // Tetgen uses 1-based enumeration
 Top_boundary_id=3;
 box_facet_boundary_id[2] = Top_boundary_id+1; 
 box_facet_boundary_id[3] = 5;
 box_facet_boundary_id[4] = 6;
 // Tetgen uses 1-based enumeration
 Bottom_boundary_id=6;
 box_facet_boundary_id[5] = Bottom_boundary_id+1;

 //Make the outer boundary object
 Outer_boundary_pt = 
  new TetgenMeshFacetedSurface(box_point,box_facet,box_facet_boundary_id);
 

 // Disk
 //-----
 unsigned npts_sheet=12;
 double sheet_radius=1.0;
 Vector<Vector<double> > sheet_point(npts_sheet); 
 for (unsigned j=0;j<npts_sheet;j++)
  {
   double phi=2.0*MathematicalConstants::Pi*double(j)/double(npts_sheet);
   sheet_point[j].resize(3);
   sheet_point[j][0]=sheet_radius*cos(phi);
   sheet_point[j][1]=sheet_radius*sin(phi);
   sheet_point[j][2]=0.0;
  }    

 //Set up the connectivity
 Vector<Vector<unsigned> > sheet_facet(1);
 sheet_facet[0].resize(npts_sheet);
 for (unsigned j=0;j<npts_sheet;j++)
  {
   sheet_facet[0][j] = j;
  }

 // Shet boundary ID (tetgen uses 1-based enumeration)
 Sheet_boundary_id=0;
 Vector<unsigned> sheet_facet_boundary_id(1,Sheet_boundary_id+1);
 
 //Create the inner boundary object
 Inner_boundary_pt.resize(1);
 Inner_boundary_pt[0] = new TetgenMeshFacetedSurface(sheet_point,
                                                     sheet_facet,
                                                     sheet_facet_boundary_id);


 // Build the mesh
 //---------------
 
 // Initial element volume
 double initial_element_volume=0.1;

 Problem::mesh_pt() = 
  new RefineableTetgenMesh<ELEMENT>(Outer_boundary_pt,
                                    Inner_boundary_pt,
                                    initial_element_volume,
                                    this->time_stepper_pt(),
                                    true);

 // Show us the mesh!
 mesh_pt()->output("orig_mesh_coarse.dat",2);  
 mesh_pt()->output("orig_mesh.dat"); 
 mesh_pt()->output_boundaries("orig_mesh_boundaries.dat");

 // Set error estimator for bulk mesh
 Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
 mesh_pt()->spatial_error_estimator_pt()=error_estimator_pt;

 // Set targets for spatial adaptivity
 mesh_pt()->max_permitted_error()=0.005; 
 mesh_pt()->min_permitted_error()=0.001; 
 mesh_pt()->max_element_size()=1.0;
 mesh_pt()->min_element_size()=0.001; 

 // Set the problem pointer
 mesh_pt()->problem_pt()=this;

 // Complete problem setup
 complete_problem_setup();
 
 // hierher use iterative solver and/or mumps
 
 // Setup equation numbering scheme
 oomph_info <<"Number of equations: " << assign_eqn_numbers() << std::endl; 
 
}



//========================================================================
/// Complete problem setup: No slip everywhere apart from top where
/// we impose parallel, axially traction free outflow
//========================================================================
template<class ELEMENT>
void DiskInContainerProblem<ELEMENT>::complete_problem_setup()
{
 // Set the boundary conditions for this problem 
 unsigned num_bound = mesh_pt()->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {

#ifdef DO_FLUID   

   // By default pin all three velocities
   unsigned final_index = 2;
   
   //Do not pin the z-velocity at top:
   if(ibound==Top_boundary_id)
    {
     final_index = 1;
    }
   
#else

   // Only single, scalar unknown in Poisson
   unsigned final_index = 0;

#endif

   // Pin whatever needs to be pinned
   unsigned num_nod= mesh_pt()->nboundary_node(ibound);
   for (unsigned inod=0;inod<num_nod;inod++)
    {
     for(unsigned i=0;i<=final_index;++i)
      {
        mesh_pt()->boundary_node_pt(ibound,inod)->pin(i);
      }
    }
  }

#ifdef DO_FLUID   
 
 // Loop over the elements to set up element-specific 
 // things that cannot be handled by constructor
 unsigned n_element = mesh_pt()->nelement();
 for(unsigned i=0;i<n_element;i++)
  {
   // Upcast from GeneralElement to the present element
   ELEMENT* el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));

   //Set the source function pointer
   el_pt->re_pt() = &Global_Parameters::Re;
  }

#endif

 // Build line visualiser
 build_line_visualiser();

 }

//========================================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void DiskInContainerProblem<ELEMENT>::doc_solution(const unsigned& nplot,
                                                DocInfo& doc_info)
{ 

 ofstream some_file;
 ofstream coarse_some_file;
 char filename[100];

 unsigned nb=mesh_pt()->nboundary();
 for (unsigned b=0;b<nb;b++)
  {
   sprintf(filename,"%s/boundary_coordinate%i_%i.dat",
           doc_info.directory().c_str(),
           b,
           doc_info.number());
   some_file.open(filename);
   mesh_pt()->Mesh::template doc_boundary_coordinates<ELEMENT>(b,some_file);
   some_file.close();
  }   

 // Output elements adjacent to disk boundary 
 sprintf(filename,"%s/elements_next_to_disk%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 sprintf(filename,"%s/coarse_elements_next_to_disk%i.dat",
         doc_info.directory().c_str(),
         doc_info.number());
 coarse_some_file.open(filename);
 unsigned n_el=mesh_pt()->nboundary_element(Sheet_boundary_id);
 oomph_info << "Number of elements next to disk: " 
            << n_el << std::endl;
 for (unsigned e=0;e<n_el;e++)
  {
   mesh_pt()->boundary_element_pt(Sheet_boundary_id,e)->
    output(some_file,nplot);
   mesh_pt()->boundary_element_pt(Sheet_boundary_id,e)->
    output(coarse_some_file,2);
  }
 some_file.close();
 coarse_some_file.close();

 // Output boundaries
 //------------------
 sprintf(filename,"%s/boundaries%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output_boundaries(some_file);
 some_file.close();


 // Output solution
 //----------------
 sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,nplot);
 some_file.close();

 // Output solution showing element outlines
 //-----------------------------------------
 sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 mesh_pt()->output(some_file,2);
 some_file.close();


 // Output line visualiser solution
 //--------------------------------
 sprintf(filename,"%s/line_soln%i.dat",doc_info.directory().c_str(),
         doc_info.number());
 some_file.open(filename);
 LV_pt->output(some_file);
 some_file.close();

} // end of doc




//========================================================================
/// Driver
//========================================================================
int main(int argc, char* argv[])
{
 // Label for output
 DocInfo doc_info;
 
 // Output directory
 doc_info.set_directory("RESLT");
  
 // Number of output points per edge
 unsigned nplot=5;


#ifdef DO_FLUID

 // Build problem 
 DiskInContainerProblem<ProjectableCrouzeixRaviartElement<
  TCrouzeixRaviartElement<3> > > problem;
 
#else

 // Build problem
 DiskInContainerProblem<ProjectablePoissonElement<
  TPoissonElement<3,3> > > problem;

#endif

 //Output initial guess
 problem.doc_solution(nplot,doc_info);
 
 //Increment counter for solutions 
 doc_info.number()++;

 // Solve the problem with one adaptation
 unsigned max_adapt=5;
 problem.steady_newton_solve(max_adapt);
 
 //Output solution
 problem.doc_solution(nplot,doc_info);
 
 //Increment counter for solutions 
 doc_info.number()++;

}



